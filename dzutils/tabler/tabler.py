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
    ):
        """
        """
        self._data_dict = {}
        self._gp_dict = gp.Dict(getpy_keytype, getpy_valuetype)
        self.df_max_rows = df_max_rows
        self.storage_dir = storage_dir
        self.data_dir = f"{storage_dir}/data_tables"
        if dataframe:
            # check the size of the dataframe
            self.add_data(dataframe)

    def add_data(self, data):
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

        for i,df in enumerate(dfs,1):
