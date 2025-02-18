from scholarly import scholarly
import requests
import re

# Buscar artículo en Google Scholar
search_query = scholarly.search_pubs('FastMLST: A Multi-core Tool for Multilocus Sequence Typing of Draft Genome Assemblies')
article = next(search_query)

citas = article['num_citations']

# Crear la URL para el badge con el número de citas
badge_url = f"https://img.shields.io/badge/citations-{citas}-blue"

# Leer el contenido del archivo README.md
with open("README.md", "r") as f:
    readme_content = f.read()

# Reemplazar el badge de citas en el archivo README.md
new_readme_content = re.sub(r'!\[Citations\]\(.*\)', f"![Citations]({badge_url})", readme_content)

# Escribir el contenido actualizado de nuevo en el archivo README.md
with open("README.md", "w") as f:
    f.write(new_readme_content)
