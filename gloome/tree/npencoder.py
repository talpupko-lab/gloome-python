import json
import numpy as np


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            obj = obj.tolist()
        if isinstance(obj, list):
            result_obj = []
            for i in obj:
                if isinstance(i, np.floating):
                    result_obj.append(float(i))
                elif isinstance(i, np.integer):
                    result_obj.append(int(i))
                else:
                    result_obj.append(i)
            return result_obj
        if isinstance(obj, list):
            result_obj = []
            for i in obj:
                if isinstance(i, np.floating):
                    result_obj.append(float(i))
                elif isinstance(i, np.integer):
                    result_obj.append(int(i))
                else:
                    result_obj.append(i)
            return tuple(result_obj)
        return super().default(obj)
