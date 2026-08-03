import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from habanero import cn
import json
import re
from . import entities as ent
from . import bakta
from . import sum_stats as st
from . import go_terms as go
from . import candidates as candidates
from itertools import islice
import requests
import time
import ast
import os
from goatools.obo_parser import GODag
from concurrent.futures import ThreadPoolExecutor
from functools import partial
import time
from time import perf_counter
from datetime import timedelta
from tqdm import tqdm

