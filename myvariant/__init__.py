'''
Python Client for MyVariant.Info services
'''
from biothings_client import get_client

__version__ = '0.3.1'

class MyVariantInfo(get_client('variant', instance=False)):
    '''This is the client for MyVariant.info web services.
    Example:

        >>> mv = MyVariantInfo()

    '''
    pass

def format_hgvs(chrom, pos, ref, alt):
    mv = MyVariantInfo()
    return mv.format_hgvs(chrom=chrom, pos=pos, ref=ref, alt=alt)
format_hgvs.__doc__ = MyVariantInfo.format_hgvs.__doc__

def get_hgvs_from_vcf(input_vcf):
    mv = MyVariantInfo()
    return mv.get_hgvs_from_vcf(input_vcf=input_vcf)
get_hgvs_from_vcf.__doc__ = MyVariantInfo.get_hgvs_from_vcf.__doc__
