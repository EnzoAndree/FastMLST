from scholarly import scholarly
import requests

# Buscar artículo en Google Scholar
search_query = scholarly.search_pubs('FastMLST: A Multi-core Tool for Multilocus Sequence Typing of Draft Genome Assemblies')
article = next(search_query)
print(article)
citas = article['num_citations']

# Crear la URL para el badge con el número de citas
badge_url = f"https://img.shields.io/badge/citations-{citas}-blue"

# Guardar el badge en el repositorio o actualizar el archivo README.md
with open("citations_badge.md", "w") as f:
    f.write(f"![Citations]({badge_url})")
