import pandas as pd
from chembl_webresource_client.new_client import new_client

target = new_client.target

target_query = target.search('Attention deficit hyperactivity disorder')
targets = pd.DataFrame.from_records(target_query)
print(targets)
