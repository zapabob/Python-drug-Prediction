import pandas as pd
from chembl_webresource_client.new_client import new_client

def extract_info(target_id):
    # Create a new client for the target
    target = new_client.target

    # Search for the target
    target_query = target.search(target_id)

    # Convert the query result into a dataframe
    targets = pd.DataFrame.from_dict(target_query)

    # Get the chembl id from the dataframe
    selected_target = targets.loc[0, 'target_chembl_id']

    # Create a new client for the activity
    activity = new_client.activity

    # Filter the activities by the selected target and 'IC50'
    res = activity.filter(target_chembl_id=selected_target).filter(standard_type='IC50')

    # Convert the result into a dataframe
    activities = pd.DataFrame.from_dict(res)

    # Keep only the rows with non-null standard_value
    df = activities[activities['standard_value'].notna()]

    # Select certain columns
    selected_columns = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    df2 = df[selected_columns]

    # Return the dataframe
    return df2

# List of target IDs
target_ids = ['CHEMBL238', 'CHEMBL239', 'CHEMBL240']

# Loop over the target IDs and call the function
for target_id in target_ids:
    df = extract_info(target_id)
    print(df)

