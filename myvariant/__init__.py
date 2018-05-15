'''
Python Client for MyVariant.Info services
'''
from biothings_client import get_client

class MyVariantInfo(get_client('variant', instance=False)):
    '''This is the client for MyVariant.info web services.
    Example:

        >>> mv = MyVariantInfo()

    '''
    pass
