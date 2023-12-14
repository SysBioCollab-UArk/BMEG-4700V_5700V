import json
import numpy as np

# Example fitted_params with potentially non-serializable objects
fitted_params = {
    "param1": 42,
    "param2": np.array([1, 2, 3]),  # NumPy array is not directly serializable
    "param3": "serializable string",
}

# Step 1: Create a function to handle serialization of non-serializable objects
def serialize_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    # Handle other non-serializable objects here if needed
    raise ValueError(f"Object of type {type(value)} is not serializable.")

# Step 2: Check and serialize the contents of fitted_params
serialized_params = {}
for key, value in fitted_params.items():
    try:
        serialized_value = serialize_value(value)
        serialized_params[key] = serialized_value
    except ValueError as e:
        print(f"Error serializing object with key '{key}': {str(e)}")

# Step 3: Save the serialized fitted_params to a JSON file
with open("fitted_params.json", "w") as json_file:
    json.dump(serialized_params, json_file, indent=4)

print("fitted_params has been saved to fitted_params.json")
