#!/bin/bash

# Verificar si se proporciona un argumento de archivo
if [ $# -ne 1 ]; then
    echo "Uso: $0 <archivo>"
    exit 1
fi

# Nombre del archivo de datos
data_file="$1"

# Verificar si el archivo de datos existe
if [ ! -f "$data_file" ]; then
    echo "El archivo $data_file no existe."
    exit 1
fi

OUTPUT_DIR=$(head -n 1 $data_file | awk '{print $1}' | sed 's|./||')
echo "OUTPUT_DIR: $OUTPUT_DIR"
# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

line_number=1
while IFS= read -r line; do
  echo "line_number: $line_number"
  # Extract file paths, ignoring the first element "./j3"
  files=($(echo "$line" | cut -d' ' -f2-))
  # Iterate over the extracted files
  file_index=1
  for file in "${files[@]}"; do
    if [[ -f "$file" ]]; then
      # Contar el número de líneas en el archivo
      line_count=$(wc -l < "$file")
      echo "file: $file : $line_count lines"

      # Save priorities to a file in the j3 directory
      output_file="$OUTPUT_DIR/pri${file_index}-${line_number}"

      # Generar una secuencia de números de 1 a line_count, luego permutarla
      #shuf -i 1-"$line_count" > "$output_file"
      # Generar una permutación de los números del 1 al line_count
      # Generar una permutación de los números del 1 al line_count y aplicar la exponencial

      shuf -i 1-"$line_count" | while read -r x; do
        div=$(echo "0.1 * $line_count" | bc -l)  # Ensure division is done correctly
        exp_value=$(awk -v x="$x" -v div="$div" 'BEGIN {print exp(x/div)}')
        #exp_value_int=$(awk -v exp_value="$exp_value" 'BEGIN {print int(exp_value + 0.5)}')
        #exp_value=$(echo "e($x/$div)" | bc -l)  # Aplicar e^x
        #exp_value=$(printf "%.0f" $(echo "$exp_value"))

        echo "$exp_value" >> "$output_file"
      done

      ((file_index++))

    else
      echo "Warning: File '$file' not found. Skipping..."
    fi
  done
  ((line_number++))
done < "$data_file"

# Contar el número de líneas en el archivo
#line_count=$(wc -l < "$data_file")
#
## Crear el archivo de prioridades
#priorities_file="test$num"
#rm -f "$priorities_file"  # Eliminar el archivo de prioridades si ya existe
#
## Generar una secuencia de números de 1 a line_count, luego permutarla
#shuf -i 1-"$line_count" > "$priorities_file"

echo "Se ha creado el archivo de prioridades \"$output_file\" con una permutación única."
