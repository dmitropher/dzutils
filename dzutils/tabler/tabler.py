import getpy as gp

# Class holding a getpy dict( with option to write to disk)
# As well as a dict with table_index:path/to/table.json
# Uses key2table and key2index


class Tabler:
    def __init__(
        self,
        getpy_keytype,
        getpy_valuetype,
        dataframe=None,
        df_max_rows=1000,
        storage_dir="./tabler_data",
        name="tably",
        gp_dict_name="dict",
    ):
        """
        """
        self._data_dict = {}
        self._gp_dict = gp.Dict(getpy_keytype, getpy_valuetype)
        self._gp_dict_name = gp_dict_name
        self.df_max_rows = df_max_rows
        self.storage_dir = storage_dir
        self.data_dir = f"{storage_dir}/data_tables"
        self.data_dir = f"{storage_dir}/dict_bin"
        self.name = name
        if dataframe:
            # check the size of the dataframe
            self.add_data(dataframe)

    def add_data(self, data, key):
        df_size = len(data.index)
        basename = data.name() if data.name() else "data"

        dfs = [data]
        if df_size > self.df_max_rows:
            """
            Chunk em, dict em, dump em
            """
            nchunks = df_size // self.df_max_rows + 1
            dfs = [
                dfs[0][i * self.df_max_rows : (i + 1) * self.df_max_rows]
                for i in range(nchunks)
            ]
        n_entries = len(self._data_dict)

        for i, df in enumerate(dfs, 1):
            """
            """
            # pick a name
            # check for duplicates
            # save to disk
            # save a reference to the file in the data dict
            # convert name to a value prefix - be mindful of the max index in the data
            # Add prefix to the index of each data point, that's values
            # Take the keys and generated values, add them to the gp dict

    def dump_dict(self, name=""):
        if not name:
            name = f"{self.name}_{self._gp_dict_name}"
        # check if the name is duplicated
        # save the gp_dict, save the name if duplicated
