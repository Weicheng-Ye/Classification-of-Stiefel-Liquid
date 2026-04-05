import os
import ast

def load_mathematica_data(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        # data.csv contains lines resembling mathematica list definitions.
        content = f.read()
    
    # Mathematica's literal lists `{1, 2}` convert nicely to python `[1, 2]`
    content = content.replace('{', '[')
    content = content.replace('}', ']')
    
    # Since it's a CSV, let's parse it as an evaluated literal
    # data.csv has a single row or multiple rows?
    # By looking at the file it was: '"[[...]]"', '"[[...]]"'
    import csv
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        parsed_rows = []
        for row in reader:
            parsed_row = []
            for element in row:
                if element.strip():
                    element_py = element.replace('{', '[').replace('}', ']')
                    parsed_row.append(ast.literal_eval(element_py))
            if parsed_row:
                parsed_rows.append(parsed_row)
        
        if len(parsed_rows) == 1:
            return parsed_rows[0]
        return parsed_rows

def main():
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(root_dir, 'data.csv')
    data_sl5_path = os.path.join(root_dir, 'dataSL5Rep.csv')

    data = load_mathematica_data(data_path)
    dataSL5Rep = load_mathematica_data(data_sl5_path)

    # Output lengths equivalently to Mathematica
    print("Number of Realizations in each case")
    
    # For[i=1, i<=4, i++] (0-index in python)
    # data[[1]][[1]][[i]] -> data[0][0][i]
    for i in range(4):
        len1 = len(data[0][0][i])
        len2 = len(data[1][0][i])
        len3 = len(data[2][0][i])
        len4 = len(dataSL5Rep[0][i])
        print(f"{len1}   {len2}   {len3}   {len4}")
        
    print()
    
    # For[i=1, i<=8, i++]
    for i in range(8):
        len1 = len(data[0][1][i])
        len2 = len(data[1][1][i])
        len3 = len(data[2][1][i])
        len4 = len(dataSL5Rep[1][i])
        print(f"{len1}   {len2}   {len3}   {len4}")

if __name__ == "__main__":
    main()
