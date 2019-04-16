"""Package Constant Values"""


VIRUSES_TAXID = 10239
LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s ' \
             '[in %(filename)s:%(lineno)d]'
CENTRIFUGE = 'centrifuge'
KRAKEN2 = 'kraken2'
classification_methods = {CENTRIFUGE, KRAKEN2}
