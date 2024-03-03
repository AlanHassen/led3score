from led3_score.fingerprints.fingerprints import Fingerprint
import pandas


class FingerprintGenerator:
    """A class to generate fingerprints
    """
    
    def __init__(self, fingerprint: Fingerprint):
        self.fingerprint = fingerprint
    
    def create_fingerprints_list(self, smiles:list):
        """
        Create a fingerprint from a smiles list
        Args:
            smiles: list of smiles molecules

        Returns:
            fingerprint
        """
        fingerprints = []
        for smile in smiles:
            fingerprints.append(self.fingerprint.create_fingerprint(smile))
        return fingerprints

    def create_fingerprints_pandas_dataframe(self, data_frame: pandas.DataFrame, smiles_column: str, fingerprint_column: str):
        """
        Create a fingerprint into a new column of the data frame
        Args:
            data_frame: data frame with smiles molecules

        Returns:
            fingerprint
        """
        fingerprints = data_frame.apply(
            lambda row: self.fingerprint.create_fingerprint(
                row[smiles_column]
            ),
            axis=1,
        )
        data_frame[fingerprint_column] = fingerprints
        return data_frame

    def flatten_fingerprints_in_dataframe(self, dataframe: pandas.DataFrame, fingerprint_column: str, fingerprint_prefix: str = "fp_"):
        """Flatten the numpy fingerpritn array into different columns

        Args:
            dataframe (pandas.DataFrame): the dataframe containing the fingerprint column
            fingerprint_column (str): the name of the fingerprint column

        Returns:
            pandas.DataFrame: the dataframe with the flattened fingerprint
        """

        # make pandas dataframe out of the numpy fingerprint
        fingerprints = dataframe.apply(lambda row: pandas.Series(row[fingerprint_column]), axis = 1)
        # rename the columns to fp + index
        fingerprints.columns = [fingerprint_prefix + str(i) for i in range(fingerprints.shape[1])]

        flattened_dataframe = pandas.concat([dataframe, fingerprints], axis = 1)

        self._assert_equal_fingerprints(flattened_dataframe)

        # drop the fingerprint column
        flattened_dataframe = flattened_dataframe.drop(columns = ["fingerprint"])

        return flattened_dataframe

    def _assert_equal_fingerprints(self, dataframe: pandas.DataFrame, fingerprint_prefix: str = "fp_"):
        """Assert that the fingerprints are equal"""
        
        fingerprint_column_df = self.get_fingerprint_columns(dataframe, fingerprint_prefix)
        
        # assert that the original fingerprint column and the flattened fingerprint columns are the same
        for index, fingerprint_row in fingerprint_column_df.iterrows():
            assert (fingerprint_row.to_numpy() == dataframe["fingerprint"].iloc[index]).all()

        print("Fingerprints are equal")

    @staticmethod
    def get_fingerprint_columns(dataframe: pandas.DataFrame, fingerprint_prefix: str = "fp_"):
        """Find the fingerprint columns in the dataframe

        Args:
            dataframe (pandas.DataFrame): the dataframe to search

        Returns:
            list: the list of fingerprint columns
        """
         # get the columns with fp_ prefix
        fingerprint_columns = [column for column in dataframe.columns if column.startswith(fingerprint_prefix)]
        
        # select the columns with fp_ prefix
        fingerprint_column_df = dataframe[fingerprint_columns]
        return fingerprint_column_df