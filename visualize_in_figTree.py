#!/usr/bin/env python3
# coding: utf-8

#################################################################################
#                                                                               #
# Script to retrieve and add taxonomic annotation to newick trees in a          #
#     FigTree-happy nexus format!                                               #
#                                                                               #
#                                         Thiberio Rangel, lthiberiol@gmail.com #
#                                                                               #
#################################################################################

import ete3
import re
import os
import sys
import argparse
import pandas as pd
import requests

#
# parse user provided inputs
#
parser = argparse.ArgumentParser(description='Add taxonomic tags to nexus leaves for a happy FigTree visualization!')

group     = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-n',
                   action='store_true',
                   help  ='Use scientific name from taxa at the start of leaf names\n'
                          '<scientific_name_whatever_comes_afterwards>')
group.add_argument('-p',
                   action='store_true',
                   help  ='Use protein accession number at the start of leaf names\n'
                          '<proteinAcc_whatever_comes_afterwards>')
group.add_argument('-e',
                   action='store_true',
                   help  ="Leaves named in eggNOG's style\n"
                          '<taxID.locus_tag>')

group2     = parser.add_argument_group('Annotation features options',
                                       'Annotated features must be provided as tables '
                                       'where the first column contain leaf names, exactly as in the tree, and '
                                       'the first row contain column names. '
                                       'Accepted formats are CSV (-c), TSV (-t), and Excel (-e).')

feat_group = group2.add_mutually_exclusive_group(required=False, )
feat_group.add_argument('-c',
                        type=str,
                        help="CSV file containing features to be added to resulting FigTree's nexus")
feat_group.add_argument('-t',
                        type=str,
                        help="Tab-separated file containing features to be added to resulting FigTree's nexus")
feat_group.add_argument('-w',
                        type=str,
                        help="Excel file containing features to be added to resulting FigTree's nexus")

parser.add_argument('-r', 
                    required=False,
                    type    =str, 
                    help    ='Path to newick containing reference rooting position to match.')

parser.add_argument('newick_file', 
                    type=str, 
                    help='newick file to be visualize in figtree')

arguments      = parser.parse_args()
by_sci_name    = arguments.n
by_protein_acc = arguments.p
as_eggNOG      = arguments.e
filename       = arguments.newick_file

if   arguments.c:
    feature_df = pd.read_csv(arguments.c, sep=',', index_col=0)
elif arguments.t:
    feature_df = pd.read_csv(arguments.t, sep='\t', index_col=0)
elif arguments.w:
    feature_df = pd.read_excel(arguments.w, index_col=0)
else:
    feature_df = pd.DataFrame()

#################################################################################

ncbi         = ete3.NCBITaxa()

if os.path.isfile(filename):
    tree = ete3.Tree(filename, format=1)
else:
    raise SystemExit('*ERROR: provided path to newick file does not exist!')

#################################################################################

if arguments.r:
    
    if os.path.isfile(arguments.r):
        ref_tree = ete3.Tree(arguments.r, format=1)
    else:
        raise SystemExit('*ERROR: provided path to newick file containing reference root does not exist!')
    
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            leaf = tree.get_leaves_by_name(node.name)[0]
            tree.set_outgroup(leaf)

        else:
            is_it_monophyletic, clade_type, fucking_up = tree.check_monophyly(
                node.get_leaf_names(), 
                'name',
                unrooted=False
            )
            if is_it_monophyletic:
                equivalent = tree.get_common_ancestor(node.get_leaf_names())
                tree.set_outgroup(equivalent)
            else:
                tree.set_outgroup(fucking_up.pop())
                equivalent = tree.get_common_ancestor(node.get_leaf_names())
                tree.set_outgroup(equivalent)
        
#################################################################################

#
# leaves are names as eggNOG
#
if as_eggNOG:
    acc2tax    = {}
    lineage_df = pd.DataFrame()
    for leaf in tree.get_leaf_names():
        tax_id        = int(leaf.split('.')[0])
        acc2tax[leaf] = tax_id

        tmp_lineage = pd.Series({rank : taxon
                                 for taxon, rank in ncbi.get_rank(
                                     ncbi.get_lineage(tax_id)).items()
                                })
        tmp_lineage = pd.Series(index=tmp_lineage.index,
                                data =ncbi.translate_to_names(tmp_lineage))

        tmp_lineage.name = tax_id

        lineage_df = lineage_df.append(tmp_lineage)

    lineage_df.drop(columns=['no rank'], inplace=True)
    lineage_df = lineage_df[~lineage_df.index.duplicated()]

#
# if starting with species name
#
if by_sci_name:
    # genera     = []
    acc2tax    = {}
    lineage_df = pd.DataFrame()
    for leaf in tree.get_leaf_names():
        regex = re.search('^(?:Candidatus_)*([A-Z][a-z]+(?:_[a-z]+\.?)?)',
                        leaf)
        if not regex:
            continue

        tmp_species     = regex.group(1).replace('_', ' ')

        translated_name = ncbi.get_name_translator([tmp_species])
        if not translated_name and tmp_species.endswith(' sp.'):
            translated_name = ncbi.get_name_translator([tmp_species.replace(' sp.', '')])

        if not translated_name:
            continue

        tax_id               = list(translated_name.values())[0][0]
        acc2tax[tmp_species] = tax_id

        tmp_lineage = pd.Series({rank : taxon
                                 for taxon, rank in ncbi.get_rank(
                                     ncbi.get_lineage(tax_id)).items()
                                })
        tmp_lineage = pd.Series(index=tmp_lineage.index,
                                data =ncbi.translate_to_names(tmp_lineage))

        tmp_lineage.name = tax_id

        lineage_df = lineage_df.append(tmp_lineage)

    lineage_df.drop(columns=['no rank'], inplace=True)
    lineage_df = lineage_df[~lineage_df.index.duplicated()]

#
# if starting with protein accession
#
if by_protein_acc:

    protein_acc = [re.match('((?:\w{2,3}_)?[^_.]+)',
                            leaf).group(1)
                   for leaf in tree.get_leaf_names()]

    acc2tax    = {}
    lineage_df = pd.DataFrame()
    for window_start in range(0, len(protein_acc), 100):
        tmp_accessions = protein_acc[window_start:window_start+100]
        annotations    = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
                                      'db=protein&'
                                      'retmode=text&'
                                      'rettype=gp&'
                                      f'id={",".join(tmp_accessions)}').text

        for block in annotations.split('\n//\n'):
            block = block.strip()

            if not block:
                continue

            tmp_acc, tmp_tax = re.search('^LOCUS\s+(\S+).+db_xref="taxon:(\d+)"',
                                         block,
                                         re.S).groups()

            acc2tax[tmp_acc] = int(tmp_tax)

            tmp_lineage = pd.Series({rank : taxon
                                     for taxon, rank in ncbi.get_rank(
                                         ncbi.get_lineage(int(tmp_tax))).items()
                                    })
            tmp_lineage = pd.Series(index=tmp_lineage.index,
                                    data =ncbi.translate_to_names(tmp_lineage))

            tmp_lineage.name = int(tmp_tax)

            lineage_df = lineage_df.append(tmp_lineage)

    lineage_df.drop(columns='no rank', inplace=True)
    lineage_df = lineage_df[~lineage_df.index.duplicated()]

#
# build FigTree happy nexus file!
#
out  = open(f'{filename}.figTree', 'w')
out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(tree))

branch_names = {}
for count, node in enumerate(tree.traverse()):
    
    if node.is_leaf():
        if by_protein_acc:
            tmp_acc = re.match('((?:\w{2,3}_)?[^_.]+)', node.name).group(1)
        elif by_sci_name:
            regex = re.search('^(?:Candidatus_)*([A-Z][a-z]+(?:_[a-z]+\.?)?)', 
                              node.name)
            tmp_acc = None
            if regex:
                tmp_acc = regex.group(1).replace('_', ' ')
        elif as_eggNOG:
            tmp_acc = node.name
        else:
            continue

        if tmp_acc not in acc2tax or            acc2tax[tmp_acc] not in lineage_df.index:
            out.write(f'\t{node.name}\n')
            continue

        tmp_tax_id = acc2tax[tmp_acc]
        
        out.write(f'\t{node.name} [&')
        for rank, taxon in lineage_df.loc[tmp_tax_id].items():
            rank = rank.replace(' ', '_')
            out.write(f'tax_{rank}="{taxon}" ')

        if node.name in feature_df.index:
            for column in feature_df.columns:
                out.write(f'feat_{column}="{feature_df.loc[node.name, column]}" ')

        out.write(']\n')

    else:
        if '/' in node.name:
            support_values = node.name.split('/')
            branch_names[f'_branch_{count}_'] = f'&1st_support={support_values[0]},2nd_support={support_values[1]}'
        else:
            branch_names[f'_branch_{count}_'] = f'&support={node.name}'

        node.name = '_branch_%i_' % count


newick_text = tree.write(format=1, dist_formatter='%.10f')
for key, value in branch_names.items():
    newick_text = newick_text.replace(key, '[%s]' % value)
out.write(';\nend;\n')
out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
out.close()

