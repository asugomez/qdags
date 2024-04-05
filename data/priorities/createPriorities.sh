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

# Contar el número de líneas en el archivo
line_count=$(wc -l < "$data_file")

# Crear el archivo de prioridades
priorities_file="priorities"
rm -f "$priorities_file"  # Eliminar el archivo de prioridades si ya existe

# Generar números aleatorios y escribirlos en el archivo de prioridades
for ((i = 0; i < line_count; i++)); do
    random_priority=$((RANDOM % line_count))  # Generar un número aleatorio entre 0 y el número de líneas
    echo "$random_priority" >> "$priorities_file"
done

echo "Se ha creado el archivo de prioridades \"$priorities_file\"."
