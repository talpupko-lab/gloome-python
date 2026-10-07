import json
import numpy as np


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)

        if isinstance(obj, np.ndarray):
            return obj.tolist()

        if isinstance(obj, (list, tuple)):
            try:
                arr = np.asarray(obj)
                if np.issubdtype(arr.dtype, np.number):
                    py_list = arr.tolist()
                    return tuple(py_list) if isinstance(obj, tuple) else py_list
            except Exception:
                pass

            res = [float(i) if isinstance(i, np.floating) else (int(i) if isinstance(i, np.integer) else i)
                   for i in obj]
            return tuple(res) if isinstance(obj, tuple) else res

        return super().default(obj)
