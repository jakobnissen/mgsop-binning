# Enforce explicit decision of output directory
if '--directory' not in sys.argv and '-d' not in sys.argv and 'binning_workdir' not in config:
    print(sys.argv)
    raise KeyError('Binning working directory not specified.')

# Parse key chain_modules
try:
    CHAIN_MODULES = bool(config['chain_modules'])
except KeyError:
    print('"chain_modules" not specified in config file.')
    raise
if CHAIN_MODULES is not True and CHAIN_MODULES is not False:
    print('Configuration key "chain_modules" must be True or False')
    raise ValueError('Configuration key "chain_modules" must be True or False')
config['chain_modules'] = CHAIN_MODULES

# Parse key configfile
if config['chain_modules']:
    try:
        CONFIGFILE = config['configfile']
    except KeyError:
        print('"chain_module set to True, but "configfile" not specified in config file.')
        raise
    if not os.path.isfile(CONFIGFILE):
        print('The file pointed to by key "configfile" cannot be found.')
        raise FileNotFoundError(CONFIGFILE)

# Get information on cores per node
try:
    CORES = int(config['cores'])
except KeyError:
    print('Cores not specified in config file.')
    raise
except TypeError:
    print('Cores cannot be interpreted as integer')
    raise
if CORES < 1:
    raise ValueError('Cores must be a positive integer')
config['cores'] = CORES

# Parse datalines
samples, experiments, runs = dict(), dict(), dict()
for line in config['data']:
    fields = line.split()

    if len(fields) not in (4, 5):
        raise ValueError('Data field {} needs 4 or 5 fields.'.format(line))

    run, experiment, sample, *paths = fields

    # No cross-category non-uniqueness in names
    if run in experiments or run in samples or run in runs:
        raise ValueError('Run {} is not unique in data.'.format(run))
    if experiment in runs or experiment in samples:
        raise ValueError('Experiment {} is also a run/sample'.format(experiment))
    if sample in runs or sample in experiments:
        raise ValueError('Sample {} is also a run/experiment'.format(sample))

    # File exists
    for path in paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)


    # Add them to dictionaries
    if sample in samples:
        samples[sample].append(experiment)
    else:
        samples[sample] = [experiment]
    if experiment in experiments:
        experiments[experiment].append(run)
    else:
        experiments[experiment] = [run]
    runs[run] = paths

# Not both SE and PE files in same experiment
for experiment, runlist in experiments.items():
    if len({len(runs[run]) for run in runlist}) != 1:
        raise ValueError('Experiment {} has both SE and PE runs.'.format(experiment))

# No experiment is in two samples (runs are ensured unique above)
if len(set().union(*samples.values())) != sum(len(v) for v in samples.values()):
    raise ValueError('No experiments can be in two different samples.')

# Create IS_PE for easy lookup which experiments are paired end
IS_PE = {exp: len(runs[runlist[0]]) == 2 for exp, runlist in experiments.items()}

# Finally add them to the config
config['samples'] = samples
config['experiments'] = experiments
config['runs'] = runs

# Make sure the local scripts are present
files = ['contigAssembler.pl', 'ChildManager.pm', 'filtercontigs.py']
files = [srcdir('../scripts/' + file) for file in files]

for file in files:
    if not os.path.isfile(file):
        raise FileNotFoundError('File not found: {}'.format(file))

# Get information on whether or not to assemble contigs
try:
    ASSEMBLE_CONTIGS = config['assemble_contigs']
except KeyError:
    print('"assemble_contigs" not found in config file')
    raise
if ASSEMBLE_CONTIGS not in (True, False):
    print('assemble_contigs must be "true" or "false"')
    raise

# Get information on platform
try:
    PLATFORM = config['platform']
except KeyError:
    print('Platform information cannot be found in config file')
    raise
if len(PLATFORM) > 70:
    print('Platform information is too long (keep below 70 characters)')
    raise ValueError('Platform information is too long (keep below 70 characters)')

# Get information on minimum gene length
try:
    MIN_CONTIG_LENGTH = int(config['min_contig_length'])
except KeyError:
    print('Minimum contig length not specified in config file')
    raise
except ValueError:
    print('Minimum contig length cannot be interpreted as integer.')
    raise

# Get information on whether to remove duplicates
try:
    REMOVE_DUPLICATES = config['remove_duplicates']
except KeyError:
    print('remove_duplicates not specified in config file.')
    raise
if REMOVE_DUPLICATES is not True and REMOVE_DUPLICATES is not False:
    raise ValueError('remove_duplicates must be "true" or "false" in config file.')

# Get information on programs to be used
names = ['BWA', 'Samtools', 'JGI', 'Metabat']
keys = tuple(name.lower() + '_path' for name in names)

for name, key in zip(names, keys):
    if key not in config:
        raise KeyError('{} path not specified in config file'.format(name))

    path = config[key]
    if not os.path.isfile(path):
        raise FileNotFoundError('{} is not an existing file'.format(path))

print('Parsed binning config file.')
