def read_parameters(file_path):
    params = {}

    with open(file_path, 'r') as f:
        for line in f:
            # Strip any leading/trailing whitespace
            line = line.strip()
            
            # Split by '='
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()

                # Handle lists (e.g., zplus)
                if value.startswith('[') and value.endswith(']'):
                    value = value.strip('[]').split(',')
                    value = [int(v.strip()) for v in value]  # Convert list items to integers
                
                # Handle numeric values
                elif value.replace('.', '', 1).isdigit():
                    if '.' in value:
                        value = float(value)  # Convert to float if it has a decimal point
                    else:
                        value = int(value)  # Convert to int if it's a whole number
                
                # Add to the dictionary
                params[key] = value
            else:
                # Handle any cases where the line might be invalid or improperly formatted
                print(f"Skipping invalid line: {line}")

    # Create variables with the same names as keys in params
    globals().update(params)

    return params
