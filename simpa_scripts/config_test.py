from simpa.utils import PathManager

if __name__ == "__main__":
    pm = PathManager(environment_path="./path_config.env")
    print("MCX path:", pm.get_mcx_binary_path())
    print("MATLAB path:", pm.get_matlab_binary_path())
    print("Save directory:", pm.get_hdf5_file_save_path())
