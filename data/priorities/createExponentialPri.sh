#!/bin/bash

# Verificar si se proporciona un argumento de archivo
if [ $# -ne 1 ]; then
    echo "Uso: $0 <archivo>"
    exit 1
fi

# Nombre del archivo de datos
data_file="$1"
num=$2

# Verificar si el archivo de datos existe
if [ ! -f "$data_file" ]; then
    echo "El archivo $data_file no existe."
    exit 1
fi

# Contar el número de líneas en el archivo
line_count=$(wc -l < "$data_file")
echo "line_count: $line_count"

# Crear el archivo de prioridades
priorities_file="pri$num"
rm -f "$priorities_file"  # Eliminar el archivo de prioridades si ya existe

# Generar valores de prioridad usando una distribución exponencial
for ((i = 0; i < line_count; i++)); do
    priority=$(echo "e($i / $line_count)" | bc -l)  # Exponential function
    echo "$priority" >> "$priorities_file"
done

echo "Se ha creado el archivo de prioridades \"$priorities_file\"."
