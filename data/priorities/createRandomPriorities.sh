#!/bin/bash
#script_dir=$(dirname "$BASH_SOURCE")
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Verificar si se proporciona un argumento de archivo
# Verificar si se proporciona un nombre de archivo
if [ $# -ne 2 ]; then
    echo "Sin archivo: $0"
    exit 1
fi

# Nombre del archivo de datos
data_file="$1"

name="$2"


# Verificar si el archivo de datos existe
if [ ! -f "$data_file" ]; then
    echo "El archivo $data_file no existe."
    exit 1
fi

# Obtener la longitud del archivo
file_length=$(wc -l < "$data_file")

# Generar una semilla aleatoria
seed=$(od -An -N2 -i /dev/random | awk '{print $1}')

# Crear el archivo de prioridades
priorities_file=("$script_dir/$name")
rm -f "$priorities_file"  # Eliminar el archivo de prioridades si ya existe

# Establecer la semilla para el generador de números aleatorios
RANDOM=$seed

# Generar números aleatorios y escribirlos en el archivo de prioridades
for ((i = 0; i < file_length; i++)); do
    random_priority=$((RANDOM % file_length + seed))  # Generar un número aleatorio dentro del rango [seed, seed + largo del archivo - 1]
    echo "$random_priority" >> "$priorities_file"
done

chmod a+x "$priorities_file"

echo "Se ha creado el archivo de prioridades \"$priorities_file\" con semilla aleatoria $seed."
