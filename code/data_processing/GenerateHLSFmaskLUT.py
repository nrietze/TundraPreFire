import pandas as pd

def generate_bit_combinations_table():
    mapping = {
        "Aerosol_Level": {
            "11": "High Aerosol",
            "10": "Medium Aerosol",
            "01": "Low Aerosol",
            "00": None
        },
        "Water": {
            "1": "Water",
            "0": None
        },
        "Snow_Ice": {
            "1": "Snow/Ice",
            "0": None
        },
        "Cloud_Shadow": {
            "1": "Cloud Shadow",
            "0": None
        },
        "Cloud_Shadow_Adjacent": {
            "1": "Adjacent Cloud Shadow",
            "0": None
        },
        "Cloud": {
            "1": "Cloud",
            "0": None
        },
        "Cirrus": {
            "1": "Reserved",
            "0": None
        }
    }
    
    def decode_quality_flags(value):
        binary = f"{value:08b}"  # Convert to 8-bit binary
        bits = {
            "Aerosol_Level": binary[:2],  # Bits 6-7
            "Water": binary[2],             # Bit 5
            "Snow_Ice": binary[3],          # Bit 4
            "Cloud_Shadow": binary[4],      # Bit 3
            "Cloud_Shadow_Adjacent": binary[5],  # Bit 2
            "Cloud": binary[6],             # Bit 1
            "Cirrus": binary[7]             # Bit 0
        }
        return bits
    
    def map_quality_flags(bits):
        labels = [mapping[key].get(bits[key], None) for key in bits]
        return ", ".join(filter(None, labels))
    
    # Generate the table
    rows = []
    for value in range(256):  # All possible values (0–255)
        bits = decode_quality_flags(value)
        label = map_quality_flags(bits)
        rows.append({"value": value, "label": label})
    
    return pd.DataFrame(rows)

# Generate and display the table
table = generate_bit_combinations_table()

# Save to a CSV file if needed
table.to_csv("data/tables/quality_flags_table.csv", index=False)

