#
# This file is autogenerated during plugin quickstart and overwritten during
# plugin makedist. DO NOT CHANGE IT if you plan to use plugin makedist to update 
# the distribution.
#

from setuptools import setup, find_packages

kwargs = {'author': 'Tristan A. Hearn',
 'author_email': 'tristan.a.hearn@nasa.gov',
 'classifiers': ['Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering'],
 'description': 'OpenMDAO implementation of the CADRE CubeSat design problem',
 'download_url': '',
 'entry_points': u'[openmdao.component]\nCADRE.solar.Solar_ExposedArea=CADRE.solar:Solar_ExposedArea\nCADRE.comm.Comm_GainPattern=CADRE.comm:Comm_GainPattern\nCADRE.orbit.Orbit_Dynamics=CADRE.orbit:Orbit_Dynamics\nCADRE.attitude.Attitude_Attitude=CADRE.attitude:Attitude_Attitude\nCADRE.comm.Comm_VectorECI=CADRE.comm:Comm_VectorECI\nCADRE.attitude.Attitude_Angular=CADRE.attitude:Attitude_Angular\nCADRE.reactionwheel.ReactionWheel_Power=CADRE.reactionwheel:ReactionWheel_Power\nCADRE.comm.Comm_BitRate=CADRE.comm:Comm_BitRate\nCADRE.orbit.Orbit_Initial=CADRE.orbit:Orbit_Initial\nCADRE.reactionwheel.ReactionWheel_Torque=CADRE.reactionwheel:ReactionWheel_Torque\nCADRE.reactionwheel.ReactionWheel_Dynamics=CADRE.reactionwheel:ReactionWheel_Dynamics\nCADRE.sun.Sun_PositionBody=CADRE.sun:Sun_PositionBody\nCADRE.comm.Comm_EarthsSpinMtx=CADRE.comm:Comm_EarthsSpinMtx\nCADRE.comm.Comm_DataDownloaded=CADRE.comm:Comm_DataDownloaded\nCADRE.CADRE_launch.Uniformity=CADRE.CADRE_launch:Uniformity\nCADRE.comm.Comm_GSposEarth=CADRE.comm:Comm_GSposEarth\nCADRE.comm.Comm_VectorBody=CADRE.comm:Comm_VectorBody\nCADRE.sun.Sun_PositionSpherical=CADRE.sun:Sun_PositionSpherical\nCADRE.attitude.Attitude_RotationMtx=CADRE.attitude:Attitude_RotationMtx\nCADRE.attitude.Attitude_Roll=CADRE.attitude:Attitude_Roll\nCADRE.power.Power_SolarPower=CADRE.power:Power_SolarPower\nCADRE.battery.BatteryConstraints=CADRE.battery:BatteryConstraints\nCADRE.attitude.Attitude_Torque=CADRE.attitude:Attitude_Torque\nCADRE.sun.Sun_PositionECI=CADRE.sun:Sun_PositionECI\nCADRE.power.Power_CellVoltage=CADRE.power:Power_CellVoltage\nCADRE.attitude.Attitude_RotationMtxRates=CADRE.attitude:Attitude_RotationMtxRates\nCADRE.battery.BatteryPower=CADRE.battery:BatteryPower\nCADRE.rk4.RK4=CADRE.rk4:RK4\nCADRE.CADRE_launch.GroundLOC=CADRE.CADRE_launch:GroundLOC\nCADRE.comm.Comm_AntRotationMtx=CADRE.comm:Comm_AntRotationMtx\nCADRE.power.Power_Total=CADRE.power:Power_Total\nCADRE.test.test_rk_deriv.RKTest=CADRE.test.test_rk_deriv:RKTest\nCADRE.CADRE_assembly.CADRE=CADRE.CADRE_assembly:CADRE\nCADRE.comm.Comm_EarthsSpin=CADRE.comm:Comm_EarthsSpin\nCADRE.comm.Comm_AntRotation=CADRE.comm:Comm_AntRotation\nCADRE.comm.Comm_VectorAnt=CADRE.comm:Comm_VectorAnt\nCADRE.comm.Comm_LOS=CADRE.comm:Comm_LOS\nCADRE.KS.KSComp=CADRE.KS:KSComp\nCADRE.reactionwheel.ReactionWheel_Motor=CADRE.reactionwheel:ReactionWheel_Motor\nCADRE.attitude.Attitude_Sideslip=CADRE.attitude:Attitude_Sideslip\nCADRE.parameters.BsplineParameters=CADRE.parameters:BsplineParameters\nCADRE.comm.Comm_VectorSpherical=CADRE.comm:Comm_VectorSpherical\nCADRE.battery.BatterySOC=CADRE.battery:BatterySOC\nCADRE.attitude.Attitude_AngularRates=CADRE.attitude:Attitude_AngularRates\nCADRE.comm.Comm_GSposECI=CADRE.comm:Comm_GSposECI\nCADRE.sun.Sun_LOS=CADRE.sun:Sun_LOS\nCADRE.CADRE_mdp.CADRE_Optimization=CADRE.CADRE_mdp:CADRE_Optimization\nCADRE.thermal_temperature.ThermalTemperature=CADRE.thermal_temperature:ThermalTemperature\nCADRE.CADRE_launch.CADRE_Launch=CADRE.CADRE_launch:CADRE_Launch\nCADRE.comm.Comm_Distance=CADRE.comm:Comm_Distance\n\n[openmdao.container]\nCADRE.solar.Solar_ExposedArea=CADRE.solar:Solar_ExposedArea\nCADRE.comm.Comm_GainPattern=CADRE.comm:Comm_GainPattern\nCADRE.orbit.Orbit_Dynamics=CADRE.orbit:Orbit_Dynamics\nCADRE.attitude.Attitude_Attitude=CADRE.attitude:Attitude_Attitude\nCADRE.comm.Comm_VectorECI=CADRE.comm:Comm_VectorECI\nCADRE.attitude.Attitude_Angular=CADRE.attitude:Attitude_Angular\nCADRE.reactionwheel.ReactionWheel_Power=CADRE.reactionwheel:ReactionWheel_Power\nCADRE.comm.Comm_BitRate=CADRE.comm:Comm_BitRate\nCADRE.orbit.Orbit_Initial=CADRE.orbit:Orbit_Initial\nCADRE.reactionwheel.ReactionWheel_Torque=CADRE.reactionwheel:ReactionWheel_Torque\nCADRE.reactionwheel.ReactionWheel_Dynamics=CADRE.reactionwheel:ReactionWheel_Dynamics\nCADRE.sun.Sun_PositionBody=CADRE.sun:Sun_PositionBody\nCADRE.comm.Comm_EarthsSpinMtx=CADRE.comm:Comm_EarthsSpinMtx\nCADRE.comm.Comm_DataDownloaded=CADRE.comm:Comm_DataDownloaded\nCADRE.CADRE_launch.Uniformity=CADRE.CADRE_launch:Uniformity\nCADRE.comm.Comm_GSposEarth=CADRE.comm:Comm_GSposEarth\nCADRE.comm.Comm_VectorBody=CADRE.comm:Comm_VectorBody\nCADRE.sun.Sun_PositionSpherical=CADRE.sun:Sun_PositionSpherical\nCADRE.attitude.Attitude_RotationMtx=CADRE.attitude:Attitude_RotationMtx\nCADRE.attitude.Attitude_Roll=CADRE.attitude:Attitude_Roll\nCADRE.power.Power_SolarPower=CADRE.power:Power_SolarPower\nCADRE.battery.BatteryConstraints=CADRE.battery:BatteryConstraints\nCADRE.attitude.Attitude_Torque=CADRE.attitude:Attitude_Torque\nCADRE.sun.Sun_PositionECI=CADRE.sun:Sun_PositionECI\nCADRE.power.Power_CellVoltage=CADRE.power:Power_CellVoltage\nCADRE.attitude.Attitude_RotationMtxRates=CADRE.attitude:Attitude_RotationMtxRates\nCADRE.battery.BatteryPower=CADRE.battery:BatteryPower\nCADRE.rk4.RK4=CADRE.rk4:RK4\nCADRE.CADRE_launch.GroundLOC=CADRE.CADRE_launch:GroundLOC\nCADRE.comm.Comm_AntRotationMtx=CADRE.comm:Comm_AntRotationMtx\nCADRE.power.Power_Total=CADRE.power:Power_Total\nCADRE.test.test_rk_deriv.RKTest=CADRE.test.test_rk_deriv:RKTest\nCADRE.CADRE_assembly.CADRE=CADRE.CADRE_assembly:CADRE\nCADRE.comm.Comm_EarthsSpin=CADRE.comm:Comm_EarthsSpin\nCADRE.comm.Comm_AntRotation=CADRE.comm:Comm_AntRotation\nCADRE.comm.Comm_VectorAnt=CADRE.comm:Comm_VectorAnt\nCADRE.comm.Comm_LOS=CADRE.comm:Comm_LOS\nCADRE.KS.KSComp=CADRE.KS:KSComp\nCADRE.reactionwheel.ReactionWheel_Motor=CADRE.reactionwheel:ReactionWheel_Motor\nCADRE.attitude.Attitude_Sideslip=CADRE.attitude:Attitude_Sideslip\nCADRE.parameters.BsplineParameters=CADRE.parameters:BsplineParameters\nCADRE.comm.Comm_VectorSpherical=CADRE.comm:Comm_VectorSpherical\nCADRE.battery.BatterySOC=CADRE.battery:BatterySOC\nCADRE.attitude.Attitude_AngularRates=CADRE.attitude:Attitude_AngularRates\nCADRE.comm.Comm_GSposECI=CADRE.comm:Comm_GSposECI\nCADRE.sun.Sun_LOS=CADRE.sun:Sun_LOS\nCADRE.CADRE_mdp.CADRE_Optimization=CADRE.CADRE_mdp:CADRE_Optimization\nCADRE.thermal_temperature.ThermalTemperature=CADRE.thermal_temperature:ThermalTemperature\nCADRE.CADRE_launch.CADRE_Launch=CADRE.CADRE_launch:CADRE_Launch\nCADRE.comm.Comm_Distance=CADRE.comm:Comm_Distance',
 'include_package_data': True,
 'install_requires': ['openmdao.main', 'MBI'],
 'keywords': ['openmdao'],
 'license': 'Apache 2.0',
 'maintainer': 'Tristan A. Hearn',
 'maintainer_email': 'tristan.a.hearn@nasa.gov',
 'name': 'CADRE',
 'package_data': {'CADRE': ['sphinx_build/html/.buildinfo',
                            'sphinx_build/html/.dummy',
                            'sphinx_build/html/full.html',
                            'sphinx_build/html/genindex.html',
                            'sphinx_build/html/glossary.html',
                            'sphinx_build/html/index.html',
                            'sphinx_build/html/launch.html',
                            'sphinx_build/html/objects.inv',
                            'sphinx_build/html/overview.html',
                            'sphinx_build/html/pkgdocs.html',
                            'sphinx_build/html/py-modindex.html',
                            'sphinx_build/html/roll.html',
                            'sphinx_build/html/search.html',
                            'sphinx_build/html/searchindex.js',
                            'sphinx_build/html/srcdocs.html',
                            'sphinx_build/html/_downloads/0_0_data.html',
                            'sphinx_build/html/_downloads/0_1_data.html',
                            'sphinx_build/html/_downloads/0_2_data.html',
                            'sphinx_build/html/_downloads/0_3_data.html',
                            'sphinx_build/html/_downloads/0_4_data.html',
                            'sphinx_build/html/_downloads/0_5_data.html',
                            'sphinx_build/html/_downloads/0_all_data.html',
                            'sphinx_build/html/_downloads/1_0_data.html',
                            'sphinx_build/html/_downloads/1_1_data.html',
                            'sphinx_build/html/_downloads/1_2_data.html',
                            'sphinx_build/html/_downloads/1_3_data.html',
                            'sphinx_build/html/_downloads/1_4_data.html',
                            'sphinx_build/html/_downloads/1_5_data.html',
                            'sphinx_build/html/_downloads/1_all_data.html',
                            'sphinx_build/html/_images/0_0.png',
                            'sphinx_build/html/_images/0_1.png',
                            'sphinx_build/html/_images/0_2.png',
                            'sphinx_build/html/_images/0_3.png',
                            'sphinx_build/html/_images/0_4.png',
                            'sphinx_build/html/_images/0_5.png',
                            'sphinx_build/html/_images/0_all.png',
                            'sphinx_build/html/_images/1_0.png',
                            'sphinx_build/html/_images/1_1.png',
                            'sphinx_build/html/_images/1_2.png',
                            'sphinx_build/html/_images/1_3.png',
                            'sphinx_build/html/_images/1_4.png',
                            'sphinx_build/html/_images/1_5.png',
                            'sphinx_build/html/_images/1_all.png',
                            'sphinx_build/html/_images/cadre3.jpg',
                            'sphinx_build/html/_images/design.png',
                            'sphinx_build/html/_images/launch.png',
                            'sphinx_build/html/_images/opt.png',
                            'sphinx_build/html/_images/roll_results.png',
                            'sphinx_build/html/_images/uniform.png',
                            'sphinx_build/html/_modules/index.html',
                            'sphinx_build/html/_modules/CADRE/attitude.html',
                            'sphinx_build/html/_modules/CADRE/battery.html',
                            'sphinx_build/html/_modules/CADRE/CADRE_assembly.html',
                            'sphinx_build/html/_modules/CADRE/CADRE_launch.html',
                            'sphinx_build/html/_modules/CADRE/CADRE_mdp.html',
                            'sphinx_build/html/_modules/CADRE/comm.html',
                            'sphinx_build/html/_modules/CADRE/kinematics.html',
                            'sphinx_build/html/_modules/CADRE/KS.html',
                            'sphinx_build/html/_modules/CADRE/orbit.html',
                            'sphinx_build/html/_modules/CADRE/parameters.html',
                            'sphinx_build/html/_modules/CADRE/power.html',
                            'sphinx_build/html/_modules/CADRE/reactionwheel.html',
                            'sphinx_build/html/_modules/CADRE/rk4.html',
                            'sphinx_build/html/_modules/CADRE/solar.html',
                            'sphinx_build/html/_modules/CADRE/sun.html',
                            'sphinx_build/html/_modules/CADRE/thermal_temperature.html',
                            'sphinx_build/html/_modules/CADRE/test/test_assembly.html',
                            'sphinx_build/html/_modules/CADRE/test/test_CADRE_derivs.html',
                            'sphinx_build/html/_modules/CADRE/test/test_components.html',
                            'sphinx_build/html/_modules/CADRE/test/test_derivatives.html',
                            'sphinx_build/html/_modules/CADRE/test/test_rk_deriv.html',
                            'sphinx_build/html/_sources/full.txt',
                            'sphinx_build/html/_sources/glossary.txt',
                            'sphinx_build/html/_sources/index.txt',
                            'sphinx_build/html/_sources/launch.txt',
                            'sphinx_build/html/_sources/overview.txt',
                            'sphinx_build/html/_sources/pkgdocs.txt',
                            'sphinx_build/html/_sources/roll.txt',
                            'sphinx_build/html/_sources/srcdocs.txt',
                            'sphinx_build/html/_static/_static',
                            'sphinx_build/html/_static/ajax-loader.gif',
                            'sphinx_build/html/_static/basic.css',
                            'sphinx_build/html/_static/comment-bright.png',
                            'sphinx_build/html/_static/comment-close.png',
                            'sphinx_build/html/_static/comment.png',
                            'sphinx_build/html/_static/default.css',
                            'sphinx_build/html/_static/doctools.js',
                            'sphinx_build/html/_static/down-pressed.png',
                            'sphinx_build/html/_static/down.png',
                            'sphinx_build/html/_static/file.png',
                            'sphinx_build/html/_static/jquery.js',
                            'sphinx_build/html/_static/minus.png',
                            'sphinx_build/html/_static/plus.png',
                            'sphinx_build/html/_static/pygments.css',
                            'sphinx_build/html/_static/searchtools.js',
                            'sphinx_build/html/_static/sidebar.js',
                            'sphinx_build/html/_static/underscore.js',
                            'sphinx_build/html/_static/up-pressed.png',
                            'sphinx_build/html/_static/up.png',
                            'sphinx_build/html/_static/websupport.js',
                            'test/__init__.py',
                            'test/data1346.pkl',
                            'test/speeds.py',
                            'test/test_assembly.py',
                            'test/test_CADRE_derivs.py',
                            'test/test_components.py',
                            'test/test_derivatives.py',
                            'test/test_rk_deriv.py']},
 'package_dir': {'': 'src'},
 'packages': ['CADRE', 'CADRE.data', 'CADRE.test'],
 'url': 'https://github.com/OpenMDAO-Plugins/CADRE',
 'version': '0.4',
 'zip_safe': False}


setup(**kwargs)

