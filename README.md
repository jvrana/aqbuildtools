# Aquarium-SD2 Strain Builder Tools (aqbt)

This repo contains many Aquarium-related tools for constructing new engineered strains.

**Features**

* Automated design and Aquarium workflow submission of plasmid assemblies
* Automated design and Aquarium workflow submission of yeast constructions
* Conversion of Aquarium strain definitions to engineered GFF with FASTA.
* JSON-schema validated serialization of Aquarium models
* Aquarium <-> Benchling integration
* Aquarium <-> SynBioHub integration (coming soon)

## Usage

### Installation

```
pip install .
```

Install BLAST locally using the following:

```
pyblast install
```

### Credentials

AqBuildTools (aqbt) requires a `credentials.toml` file to connect with 
Aquarium, Benchling and SynBioHub (coming soon). The following 
creates the `default` scope for a aqbt session:
```
[aquarium.production]
login = 'mylogin'
password = 'mypass'
url = 'http://aquarium.org'

[benchling.production]
apikey = 'sk_398dfg983nsdfsdflksj'
initials = 'UWBF'
schema = 'Aquarium DNA'
prefix = 'aq'
folder_id = 'lib_ILUNzz6N'
id = 'src_KE6uFvex'

[session.default]
aquarium = 'production'
```

```python
from aqbt import AquariumBuildTools

aqbt = AquariumBuildTools.from_toml('credentials.toml')
```

### Examples

#### Querying the Aquarium server

#### Connecting to Benchling and Synbiohub

#### Plasmid construction

#### Strain construction

##### Generating GFFs

## Related repos

[trident](https://github.com/klavinslab/trident),
[DASi](https://github.com/jvrana/DASi-DNA-Design),
[Benchling API](https://github.com/klavinslab/benchling-api),
[Terrarium](https://github.com/jvrana/Terrarium),
[pyblast](https://github.com/jvrana/pyblast),
[primer3-py-plus](https://github.com/jvrana/primer3-py-plus)
