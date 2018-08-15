"""ParameterSet

Description of the module

"""


class ParameterSet(object):
    """Class containing model paramters

    Attributes
    ----------
    pars : dict
        Supplied parameters
    defaultpars : dict
        Default parameters overwriting supplied parameters

    Methods
    -------
    update(pars)
        Adds and updates parameters.
    """


    def __init__(self, pars={}, defaultpars={}):
        mypars = pars.copy()
        for k in defaultpars.keys():
            mypars.setdefault(k, defaultpars[k])
        for k,v in mypars.items():
            setattr(self,k,v)

    def update(self, pars):
        """Adds and updates parameters

        Parameters
        ----------
        pars : dict or modelbase.parameters.ParameterSet
            Object containing new parameters

        Returns
        -------
        None

        Warns
        -----
        OverwritingKeys
            Prints warning if keys are overwritten
        """
        if isinstance(pars,dict):
            replaced_keys = [key for key in self.__dict__.keys() if key in pars]
            if replaced_keys:
                print("Warning: overwriting keys", replaced_keys)
            for k,v in pars.items():
                setattr(self,k,v)
        elif isinstance(pars,ParameterSet):
            self.update(pars.__dict__)
        else:
            raise TypeError("Function requires dict or ParameterSet input")
