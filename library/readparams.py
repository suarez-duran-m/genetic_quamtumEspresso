import json
import math

class readparams:
    def __init__(self):
        pass

    def readingparameters(self, filename="params.json", mode=1):
        with open(filename, "r") as f:
            data = json.load(f)

        # Retrieve the sub-dictionary for the given mode.
        read4mode = "singleSpice" if mode else "twoSpices"

        params_data = data.get(read4mode)
        if params_data is None:
            raise ValueError(f"Any of Mode \"singleSpice\" or \"twoSpices\" was found in {filename}.")

        # Set each parameter as an attribute.
        for key, value in params_data.items():
            setattr(self, key, value)

        # Post-process 'getZeroPopu'
        if not hasattr(self, "getZeroPopu"):
            self.getZeroPopu = None
        else:
            # If getZeroPopu is a string, check if it's empty or "nan"
            if isinstance(self.getZeroPopu, str):
                if self.getZeroPopu.strip() == "" or self.getZeroPopu.lower() == "nan":
                    self.getZeroPopu = None
            # If getZeroPopu is a number (int or float) and equals 0 or is NaN, set it to None
            elif isinstance(self.getZeroPopu, (int, float)):
                if self.getZeroPopu == 0 or (isinstance(self.getZeroPopu, float) and math.isnan(self.getZeroPopu)):
                    self.getZeroPopu = None

    def fetchparams(self, filename="params.json", mode=1):
        """Loads parameters from a JSON file using the specified mode and returns itself."""
        self.readingparameters(filename, mode)
        return self