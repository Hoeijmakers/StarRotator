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
                          'inclination','RpRs','P','veq','stelinc',
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
