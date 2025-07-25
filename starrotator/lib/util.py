def read_into_dictionary(filepath):
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            # Skip empty or whitespace-only lines
            if not line.strip():
                continue

            parts = line.strip().split(None, 1)  # Split on first whitespace
            if len(parts) != 2:
                raise ValueError(f"Line format error: {line.strip()}")
            key, val = parts

            # Try to cast to float, else keep as string
            try:
                val = float(val)
            except ValueError:
                pass

            data[key] = val
    return data


def check_integrity_input(input_dict,additional_keywords=[]):
    mandatory_keywords = ['sma_Rs','e','omega','inclination','obliquity',
                          'RpRs','P','veq','stelinc',
                          'drr','T','FeH','logg','u1','u2','mus','R','model']
    if len(additional_keywords) > 0:
        mandatory_keywords=mandatory_keywords+additional_keywords

    for keyword in mandatory_keywords:
        assert keyword in input_dict





def vartest(var,varname='',nan=False,pos=False,notnegative=False,types=[],dims=[]):
    """This is a shorthand containing the tests below to be able to test input on multiple
    criteria at once."""
    if nan:
        nantest(var,varname)
    if pos:
        postest(var,varname)
    if notnegative:
        notnegativetest(var,varname)
    if len(types) > 0:
        typetest(var,types,varname)
    if len(dims) > 0:
        dimtest(var,dims,varname)


def nantest(var,varname=''):
    import numpy as np
    if np.isnan(var).any()  == True:
        raise ValueError("Variable %s contains NaNs." % varname)
    if np.isinf(var).any()  == True:
        raise ValueError("Variable %s contains in-finite values." % varname)

def postest(a,varname=''):
    """This function tests whether a number/array is strictly positive."""
    import numpy as np
    if np.min(a) <= 0:
        raise ValueError('Variable %s is not strictly positive' % varname)

def notnegativetest(a,varname=''):
    """This function tests whether a number/array is strictly positive."""
    import numpy as np
    if np.min(a) < 0:
        raise ValueError('Variable %s is negative.' % varname)

def file_exists(file,varname=''):
    """This program tests if a file exists, and prints an error otherwise."""
    from os import path
    import sys
    typetest(file,str,varname='file in test.file_exists')
    if path.isfile(file) != True:
        print('Error: File %s does not exist.' % varname)
        sys.exit()

def dir_exists(dir,varname=''):
    """This program tests if a directory exists, and prints an error otherwise."""
    from os import path
    import sys
    typetest(dir,str,varname='dir in test.dir_exists')
    if path.isdir(dir) != True:
        print('Error: Directory %s does not exist.' % varname)
        sys.exit()

def typetest(var,vartype,varname=''):
    """This program tests the type of var which has the name varname against
    the type vartype (or the types in list vartype), and raises an exception if
    either varname is not a string or if type(var) is not equal to (any element in)
    vartype.

    Example:
    a = 'ohai'
    utils.typtest('a',a,str)"""
    if isinstance(varname,str) != True:
        print(varname)
        raise Exception("Input error in typetest: varname should be of type string.")

    #Test two cases, if vartype is a list, we test each element.
    if isinstance(vartype,list) == True:
        error = 1
        msg = "Type error: Variable %s should be one of these types: " % varname
        for type in vartype:
            msg+=' %s' % type
            if isinstance(var,type) == True:
                error = 0
        if error == 1:
            raise Exception(msg)
    else:
        #Otherwise, we do the classical typetest.
        if isinstance(var,vartype) != True:
            raise Exception("Type error: Variable %s should be  of type %s." % (varname,vartype))

def typetest_array(var,vartype,varname=''):
    """This program tests the type of the elements in the array or list var which has the
    name varname, against the type vartype, and raises an exception if either
    varname is not a string, type(var) is not equal to numpy.array or list, or the elements of
    var are not ALL of a type equal to vartype.
    """
    #NEED TO FIX: MAKE SURE THAT A RANGE OF TYPES CAN BE TESTED FOR, SUCH AS
    #float, np.float32, np.float64... should all pass as a float.
    import numpy as np
    if isinstance(varname,str) != True:
        raise Exception("Input error in typetest: varname should be of type string.")
    if (isinstance(var,list) != True) and (isinstance(var,np.ndarray) != True):
        raise Exception("Input error in typetest_array: %s should be of class list or numpy array." % varname)
    for i in range(0,len(var)):
        typetest(var[i],vartype,varname='element %s of %s' % (i,varname))

def dimtest(var,sizes,varname=''):
    """This program tests the dimensions and shape of the input array var.
    Sizes is the number of elements on each axis.
    The program uses the above type tests to make sure that the input is ok.
    If an element in sizes is set to zero, that dimension is not checked against.
    Example:
    import numpy as np
    a=[[1,2,3],[4,3,9]]
    b=np.array(a)
    dimtest(a,[2,3])
    dimtest(a,[3,10])
    """
    import numpy as np
    typetest(sizes,list,varname='sizes in dimtest')
    typetest_array(sizes,int,varname='sizes in dimtest')

    ndim=len(sizes)

    #First check that the number of dimensions is correct
    if np.ndim(var) != ndim:
        raise Exception("Dimension Error:  Variable %s ndim = %s but was required to be %s." % (varname,np.ndim(var),ndim))

    #Then check that each of the axes match.
    sizes_var=np.shape(var)
    for i in range(0,len(sizes)):
        if sizes[i] < 0:
            raise Exception("Sizes in dimtest was not set correctly. It contains negative values. (%s)" % sizes(i))
        if sizes[i] > 0:
            if sizes[i] != sizes_var[i]:
                raise Exception("Dimension Error: Axis %s of variable %s contains %s elements, but %s were required." % (i,varname,sizes_var[i],sizes[i]))


from jax import jit
import jax.numpy as jnp

@jit
def gaussian(x,A,mu,sig,cont=0.0):
    """This produces a gaussian function on the grid x with amplitude A, mean mu
    and standard deviation sig. Will need to expand it with a version that has
    a polynomial continuum in the same way that IDL does it."""
    return A * jnp.exp(-0.5*(x - mu)/sig*(x - mu)/sig)+cont








def create_cache_dir(user_defined_path = None):
    from pathlib import Path
    from platformdirs import user_cache_dir


    if user_defined_path is None:
        APP_NAME = "StarRotator"
        APP_AUTHOR = "YourName"  # or organization

        CACHE_DIR = Path(user_cache_dir(APP_NAME, APP_AUTHOR))
    else: 
        CACHE_DIR = Path(user_defined_path)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    print(f"--- Created cache directory at: {CACHE_DIR}")


# def get_cache_dir() -> Path:
#     return CACHE_DIR

# def clear_cache():
#     for file in CACHE_DIR.glob("*"):
#         file.unlink()





import json
from pathlib import Path
import platformdirs

# Harcoded paths that are mostly OS independent:
CONFIG_DIR = Path.home() / ".starrotator"
CONFIG_FILE = CONFIG_DIR / "config.json"

def get_default_cache_dir() -> Path:
    """Returns the default cache location (AppData folder or equivalent)."""
    return Path(platformdirs.user_cache_dir("StarRotator","Hoeijmakers"))


def save_default_config() -> None:
    save_config({'cache_dir' :  str(get_default_cache_dir())})

def save_config(config: dict) -> None:
    """Saves a config.json to the home directory. In the rare cases that this is not allowed,
    the non-existence of a config.json will be handled by the other cache functions in a default way."""
    typetest(config,dict,varname='config in util.save_config()')

    if "cache_dir" not in config:
        raise ValueError("config variable in util.save_config() requires a key 'cache_dir'")
    
    try: 
        CONFIG_DIR.mkdir(parents=True, exist_ok=True)
        with open(CONFIG_FILE, "w") as f:
            json.dump(config, f, indent=4)
        print(f'--- Created/overwrote a config.json file at {CONFIG_FILE}')
    except PermissionError:
        print(f"Warning: Cannot write to {CONFIG_FILE}. All cache actions will proceed with default cache location at {str(get_default_cache_dir())}.")
        print(f"To alleviate this warning and enable usage of a custom cache, provide access to {str(CONFIG_DIR)}.")


def load_config() -> dict:
    """Load config.json that points to the cache location. If it doesn't exist,
      or can't be read, a warning will be printed and the default cache 
      location will be assumed."""
    if not CONFIG_FILE.exists():
        print(f"Warning: No {CONFIG_FILE} exists. Proceeding with default cache location at {str(get_default_cache_dir())}.")
        print('To remove this warning while still using the default cache, create a default config file by running:')
        print('>>> from starrotator.lib.util import save_default_config')
        print('>>> save_default_config()')

    else:
        try:
            with open(CONFIG_FILE, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, PermissionError):
            print(f"Warning: {CONFIG_FILE} exists but cannot be read or decoded.")
            print(f"Using default cache location at {str(get_default_cache_dir())}.")
    # If we can't read or file doesn't exist, return the default cache location:
    return({"cache_dir": str(get_default_cache_dir())})


def set_cache_dir(new_path: str):
    """Update the cache directory in the config file with a user-defined path. 
    If the config file does not exist, it is created."""
    config = load_config()
    config["cache_dir"] = str(Path(new_path).expanduser())
    save_config(config)
    print(f"--- Cache directory updated to: {config['cache_dir']}")

def make_cache_dir_if_not_exists(cache_dir: Path):
        cache_dir = Path(cache_dir)
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)  # Create if missing
        except:
            print(f"Warning: Cannot create a cache directory at {str(cache_dir)}.")
            print(f"This means that functionality that relies on PHOENIX models will be unavailable.")
            print(f"To alleviate this, make sure that a cache can be created in the default location {str(get_default_cache_dir())}")
            print(f"Or provide a custom cache location using >>> util.set_cache_dir('your/custom/path/to/a/directory')")

def get_cache_dir() -> Path:
    """Get current cache directory from the config.json file"""
    config = load_config()
    cache_dir = Path(config["cache_dir"])
    make_cache_dir_if_not_exists(cache_dir)
    return cache_dir





# def clear_cache():
#     """Delete all files in the cache directory."""
#     cache_dir = get_cache_dir()
#     for file in cache_dir.glob("*"):
#         try:
#             file.unlink()
#         except IsADirectoryError:
#             for subfile in file.glob("*"):
#                 subfile.unlink()
#             file.rmdir()
#     print("[StarRotator] Cache cleared.")