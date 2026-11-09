import numpy as np
import f90nml

type_map = {
    "complex(wp)": np.complex128,
    "real(wp)": np.float64,
}


def customLoad(fname):
    descr = f90nml.read(str(fname) + "_descr.nml")["arr_descr"]
    shap = descr["shap"]
    dtype = type_map[descr["dtype"]]

    return np.fromfile(str(fname) + ".bin", dtype=dtype).reshape(shap, order="F")
