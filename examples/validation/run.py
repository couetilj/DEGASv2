"""Local runner; the same engine also supplies parallel SLURM workers and pooling."""
import sys
from DEGAS_python.sweep import main

if __name__ == '__main__':
    main(['run', *sys.argv[1:]])
