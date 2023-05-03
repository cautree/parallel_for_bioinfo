
#attributes in gff file
#'ID=gene:ENSG00000227232;Name=WASH7P;biotype=unprocessed_pseudogene;description=WAS protein family homolog 7 pseudogene [Source:HGNC Symbol%3BAcc:HGNC:38034];\
# gene_id=ENSG00000227232;\
#havana_gene=OTTHUMG00000000958;havana_version=1;logic_name=havana;version=5'

import re

RE_GENE_NAME = re.compile(r'Name=(?P<gene_name>.+?);')
def extract_gene_name(attributes_str):
    res = RE_GENE_NAME.search(attributes_str)
    return res.group('gene_name')


ndf['gene_name'] = ndf.attributes.apply(extract_gene_name)

RE_GENE_ID = re.compile(r'gene_id=(?P<gene_id>ENSG.+?);')
def extract_gene_id(attributes_str):
    res = RE_GENE_ID.search(attributes_str)
    return res.group('gene_id')


ndf['gene_id'] = ndf.attributes.apply(extract_gene_id)


RE_DESC = re.compile('description=(?P<desc>.+?);')
def extract_description(attributes_str):
    res = RE_DESC.search(attributes_str)
    if res is None:
        return ''
    else:
        return res.group('desc')


ndf['desc'] = ndf.attributes.apply(extract_description)