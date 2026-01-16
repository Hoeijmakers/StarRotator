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


def check_integrity_input(input_dict,req_keywords):
    """This tests the integrity of an input dictionary.
    The pattern is that the req_keyword dict has the required keywords
    as entries, set to a 2-element list, of which the first element
    is a list of the required variable types, and the second element
    is named tests that the keyword in the input dict should be tested against.
    
    Valid tests are 'type', 'pos', 'notnegative', 'nonans' and 'exists' (for paths)."""
    for keyword in req_keywords:
        if keyword not in input_dict:
            raise Exception(f'Keyword {keyword} not found in input dict.')
        types,tests = req_keywords[keyword]
        if 'type' in tests:
            typetest(input_dict[keyword],types,f'{keyword} in input dictionary')
        if 'notnegative' in tests:
            notnegativetest(input_dict[keyword],f'{keyword} in input dictionary')
        if 'pos' in tests:
            postest(input_dict[keyword],f'{keyword} in input dictionary')
        if 'nonans' in tests:
            nantest(input_dict[keyword],f'{keyword} in input dictionary')
        if 'exists' in tests:
            file_exists(input_dict[keyword],f'{keyword} in input dictionary')






def vartest(var,varname='',nonans=False,pos=False,notnegative=False,types=[],dims=[]):
    """This is a shorthand containing the tests below to be able to test input on multiple
    criteria at once."""
    if nonans:
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
    """This program tests if a path exists, and raises an error otherwise."""
    from pathlib import Path
    typetest(file,[str,Path],varname='file in test.file_exists')

    if Path(file).exists() != True:
        raise FileNotFoundError(f'{varname} does not exist.')
    if Path(file).is_file() != True:
        raise FileNotFoundError(f'{varname} is not a file.')

def dir_exists(dir,varname=''):
    """This program tests if a directory exists, and prints an error otherwise."""
    from pathlib import Path
    import sys
    typetest(dir,[str,Path],varname='dir in test.dir_exists')
    if Path(dir).exists() != True:
        raise FileNotFoundError(f'{varname} does not exist.')
    if Path(dir).is_dir() != True:
        raise FileNotFoundError(f'{varname} is not a directory.')

def is_array_like_non_scalar(x):
    import numpy as np
    try:
        arr = np.asarray(x)
        return arr.ndim > 0
    except Exception:
        return False

def typetest(var,vartype,varname=''):
    """This program tests the type of var which has the name varname against
    the type vartype (or the types in list vartype), and raises an exception if
    either varname is not a string or if type(var) is not equal to (any element in)
    vartype.

    Example:
    a = 'ohai'
    utils.typtest('a',a,str)"""
    import numpy as np
    if isinstance(varname,str) != True:
        raise Exception("Input error in typetest: varname should be of type string.")

    #Test two cases, if vartype is a list, we test each element.
    if isinstance(vartype,list) == True:
        error = 1
        msg = "Type error: Variable %s should be one of these types: " % varname
        for T in vartype:
            msg+=' %s' % T
            if isinstance(T,str):
                if T == 'array-like':
                    if is_array_like_non_scalar(var):
                        error = 0
                else:
                    raise Exception(f'Input error in typetest: This typestring {T} is not supported.')
            else:
                if isinstance(var,T) == True:
                    error = 0
        if error == 1:
            raise Exception(msg)
    else:
        #Otherwise, we do the classical typetest.
        if isinstance(vartype,str):
            if vartype == 'array-like':
                if not is_array_like_non_scalar(var):
                    raise Exception("Type error: Variable %s should be  of type %s." % (varname,vartype))
            else:
                raise Exception('Input error in typetest: This typestring {T} is not supported.')
        else: 
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





import json
from pathlib import Path
import platformdirs
import astropy.units as u
from astropy.constants import (G as quantity_G, c as quantity_c,  # pyright: ignore[reportAttributeAccessIssue]
                               M_sun as quantity_M_sun, # pyright: ignore[reportAttributeAccessIssue]
                               R_sun as quantity_R_sun # pyright: ignore[reportAttributeAccessIssue]
                            )#


#Hardcoded natural constants in cgs:
G = quantity_G.cgs.value
c = quantity_c.cgs.value
M_sun = quantity_M_sun.cgs.value
R_sun = quantity_R_sun.cgs.value
d_in_seconds =  (1 * u.d).cgs.value # pyright: ignore[reportAttributeAccessIssue]
rad_in_deg = (1 * u.rad).to('degree').value # pyright: ignore[reportAttributeAccessIssue]



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