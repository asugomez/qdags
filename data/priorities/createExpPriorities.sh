#!/bin/bash

# Verificar si se proporciona un argumento de archivo
if [ $# -ne 1 ]; then
    echo "Uso: $0 <archivo>"
    exit 1
fi

# Input file containing the commands
INPUT_FILE="$1" #"input.txt"

# Output directory
OUTPUT_DIR="test"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Read each line from the input file
line_number=1
while IFS= read -r line; do
    echo "line: $line"
    # Extract file paths, ignoring the first element "./j3"
    files=($(echo "$line" | cut -d' ' -f2-))
    echo "files: ${files[0]}"
    echo "files: ${files[1]}"
    echo "files: ${files[2]}"
    echo "files: ${files[3]}"
    # Iterate over the extracted files
    file_index=1

    # crear permutaión con el número de lineas del archivo
    # aplicar exp(x) a cada uno de esos números
    for file in "${files[@]}"; do
        echo "file: $file"
        if [[ -f "$file" ]]; then
            # Get the number of lines in the file
            num_lines=$(wc -l < "$file")
            echo "num_lines: $num_lines"

            # Save priorities to a file in the j3 directory
            output_file="$OUTPUT_DIR/pri${file_index}-${line_number}"

            for ((i = 0; i < num_lines; i++)); do
                priority=$(echo "e($i / $num_lines)" | bc -l)  # Exponential function
                echo "$priority" >> "$output_file"
            done

            echo "priorities: ${priorities}"

            echo "output_file: $output_file"
#            printf "%s\n" "${priorities[@]}" > "$output_file"

            ((file_index++))
        else
            echo "Warning: File '$file' not found. Skipping..."
        fi
    done
    ((line_number++))
done < "$INPUT_FILE"

echo "Priority files created in '$OUTPUT_DIR' successfully."
