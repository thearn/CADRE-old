===========
Example: Optimization of the CADRE roll angle
===========

In this example, we will optimize the roll angle of the CADRE satellite as
it passes over the ground station to maximize the gain of the communications system. For the sake of simplicity of the example, the roll angle will be the only design variable considered.

First, be sure that you have installed the 'MBI' library and CADRE plugin
for OpenMDAO. In a new python file, we then import libraries that we will use

.. code-block:: python
    from openmdao.main.api import Component
    from openmdao.main.datatypes.api import Float, Array
    from CADRE import CADRE
    from openmdao.lib.drivers.api import SLSQPdriver
    import os
    import pylab
    import numpy as np

Now, we define an OpenMDAO component that computes the sum of an input array. We will use the output of this component as our objective function
later on

.. code-block:: python

    class NetGain(Component):

        """
        Computes the sum of an inputted gain array
        """
        net = Float(iotype="out")

        def __init__(self, n):
            super(NetGain, self).__init__()
            self.n = n
            self.add('gain', Array(np.zeros(n), iotype='in', shape=(n,)))

        def execute(self):
            self.net = sum(self.gain)

        def apply_derivT(self, arg, result):
            if 'gain' in result and 'net' in arg:
                result['gain'] += arg['net'] * np.ones(self.n)

Next, we create an instance of the CADRE base assembly, and set appropriate
configuration data (starting position and velocity of the satellite, and number of days since launch)

.. code-block:: python

    n, m = 1500, 150
    top = CADRE(n, m)

    # orbit initial position and velocity
    r_e2b_I0 = [-4969.91222,  4624.84149,
                1135.9414,  0.1874654, -1.62801666,  7.4302362]

    # number of days since launch
    LD = 5417.5

    top.set("LD", LD)
    top.set("r_e2b_I0", r_e2b_I0)


Running the assembly as-is gives us a baseline state of the model, with all design variables at their default values.
Our objective is to maximize the total communication gain (as computed by the Comm_GainPattern component in the CADRE assembly),
so lets get that value:

.. code-block:: python

    # Run model to get baseline net gain value
    top.run()
    obj1 = sum(top.Comm_GainPattern.gain)
    print "Net comm gain before optimization:", obj1

Now we're ready to optimize. Replace the default "RunOnce" driver with the
`SLSQPdriver()` optimization driver, add in the NetGain component, and configure the optimization problem:

.. code-block:: python

    # Add in optimization driver
    top.add("driver", SLSQPdriver())

    top.add("NetGain", NetGain(n))
    top.driver.workflow.add("NetGain")

    top.connect("Comm_GainPattern.gain", "NetGain.gain")

    top.driver.add_parameter("CP_gamma", low=0, high=np.pi / 2.)
    top.driver.add_objective("-NetGain.net")

Run the assembly to perform the optimization, and then record the new value of the gain:

.. code-block:: python

    top.run()
    obj2 = sum(top.Comm_GainPattern.gain)

This value should be about a 23% improvement over the baseline.

We can plot the roll angle, gamma, to visualize the craft roll angle over time, as selected by the optimizer:

.. code-block:: python

    pylab.figure()
    pylab.plot(top.CP_gamma)
    pylab.show()

This is implemented in `example_roll.py`, in the top-level directory of the CADRE plugin repository.