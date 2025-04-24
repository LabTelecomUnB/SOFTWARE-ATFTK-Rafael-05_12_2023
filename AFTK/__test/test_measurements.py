#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""
import h5py
import numpy
import pytest

from __test.data import get_laurent_data

from pathlib import Path

from aftk.measurements import RadiationFieldData, rot_x, rot_y, rot_z


from matplotlib import rc

here = Path(__file__).parent


class TestRadiationFieldData:

    @pytest.fixture()
    def samples(self) -> tuple:
        antenna = 1
        freq = 17

        data = get_laurent_data(antenna, freq)

        Eh = data.Eh.values
        Ev = data.Ev.values

        theta = data.theta.values
        theta[theta == 0] = 360
        theta -= 180
        theta = numpy.deg2rad(theta)

        phi = data.phi.values - 60
        phi = numpy.deg2rad(phi)

        # ---------- ---------- ---------- ---------- ---------- ----------
        # MBA = rot_y(0*numpy.pi / 2) @ rot_x(0*numpy.pi/2)

        return theta, phi, Ev, Eh

        # return RadiationFieldData.rotate_data(MBA, theta, phi, Ev, Eh)

    def test_data_size(self, samples):

        print()

        theta, phi, E_theta, E_phi = samples

        no_pole_cond = (theta > 0) & (theta < 180)

        theta = theta[no_pole_cond]
        phi = phi[no_pole_cond]
        E_theta = E_theta[no_pole_cond]
        E_phi = E_phi[no_pole_cond]

        M = len(theta)

        l_max = int(numpy.sqrt(M + 1) - 1)
        print(f'{M = }')
        print(f'{l_max = }')

    def test_initialising_and_plotting(self, samples):

        theta, phi, E_theta, E_phi = samples
        rfd = RadiationFieldData(theta, phi, E_theta, E_phi)

        rfd.plot_measured_directions()
        rfd.plot_2D_measured_radiation_pattern(show=True)
        rfd.plot_3D_measured_radiation_pattern()

    def test_creating_and_saving_regressor_matrix(self, samples):
        print('\n')

        theta, phi, E_theta, E_phi = samples
        rfd = RadiationFieldData(theta, phi, E_theta, E_phi,
                                 name='some name',
                                 description='some description',
                                 frequency=17e9)

        rfd.create_regressors_matrix(60)
        rfd.save(here / 'testing')

    def test_estimate_mode_coeff(self, samples):
        rc('text', usetex=True)
        rc('font', family='serif')
        print('\n')

        theta, phi, E_theta, E_phi = samples
        rfd = RadiationFieldData(theta, phi, E_theta, E_phi,
                                 name='some name',
                                 description='some description',
                                 frequency=17e9)

        q, P, res, max_rel_error, \
            theta, phi, E_theta, E_phi, \
            E_est_theta, E_est_phi = rfd.estimate_mode_coeff(10)

        print(f'mean residual: {res.mean()}')
        print(max_rel_error, P.diagonal().max())

        rfd.plot_mode_power_distribution(q)

    def test_estimate_mode_coeff_unbiased_and_regularised(self, samples):
        rc('text', usetex=True)
        rc('font', family='serif')
        print('\n')

        # ========== ========== ========== ========== ========== ==========
        filepath = here / 'RULE_analysis_results.hdf5'

        if not filepath.is_file():

            theta, phi, E_theta, E_phi = samples
            rfd = RadiationFieldData(theta, phi, E_theta, E_phi,
                                     name='some name',
                                     description='some description',
                                     frequency=17e9)

            theta = rfd.theta[rfd.no_pole_cond]
            phi = rfd.phi[rfd.no_pole_cond]
            E_mea_theta = rfd.E_theta[rfd.no_pole_cond]
            E_mea_phi = rfd.E_phi[rfd.no_pole_cond]

            E = numpy.hstack([
                E_mea_theta.reshape(-1, 1),
                E_mea_phi.reshape(-1, 1)
            ]).reshape(-1, 1)

            power_rad = (E.conjugate().T @ E).real[0, 0]

            M = len(E_mea_theta)

            k_min = 10
            k_max = 70 #int(numpy.sqrt(M + 1) - 1)

            k = numpy.arange(k_min, k_max+1, 5)

            # ---------- ---------- ---------- ---------- ---------- ----------
            with h5py.File(filepath, 'a') as file:

                # ========== ========== ========== ========== ========== ==========
                measurement_group = file.create_group('measurement')

                theta_dataset = measurement_group.create_dataset('theta', data=theta)
                theta_dataset.attrs['unit'] = 'rad'
                theta_dataset.attrs['description'] = 'polar angle'

                phi_dataset = measurement_group.create_dataset('phi', data=phi)
                phi_dataset.attrs['unit'] = 'rad'
                phi_dataset.attrs['description'] = 'azimuthal angle'

                E_theta_dataset = measurement_group.create_dataset('E_theta', data=E_mea_theta)
                E_theta_dataset.attrs['unit'] = 'volt'
                E_theta_dataset.attrs['description'] = 'radiation field polar component'

                E_phi_dataset = measurement_group.create_dataset('E_phi', data=E_mea_phi)
                E_phi_dataset.attrs['unit'] = 'volt'
                E_phi_dataset.attrs['description'] = 'radiation field azimuthal component'

                # ========== ========== ========== ========== ========== ==========
                regression_group = file.create_group('regression')

                for tkp in [0, 1, 10, 100, 1000, 10000]:

                    tkp_group = regression_group.create_group(f'TKP_{tkp:05d}')

                    for _k in k:

                        k_group = tkp_group.create_group(f'k_{_k:02d}')

                        if tkp == 0:
                            q, P, res, max_rel_error, \
                                theta, phi, E_theta, E_phi, \
                                E_est_theta, E_est_phi = \
                                rfd.estimate_mode_coeff(_k)

                        else:
                            q, P, res, max_rel_error, \
                                theta, phi, E_theta, E_phi, \
                                E_est_theta, E_est_phi = \
                                rfd.estimate_mode_coeff_unbiased_and_regularised(_k, tkp)

                        residual_energy = (res.conjugate().T @ res).real[0, 0] / power_rad
                        norm_trace_cov = numpy.trace(P).real / (_k*(_k + 2))

                        k_group.attrs['residual_energy'] = residual_energy
                        k_group.attrs['norm_trace_cov'] = norm_trace_cov

                        q_dataset = k_group.create_dataset('q', data=q)
                        q_dataset.attrs['unit'] = 'watt'
                        q_dataset.attrs['description'] = 'Estimated spherical mode coeff'

                        P_dataset = k_group.create_dataset('P', data=P)
                        P_dataset.attrs['unit'] = 'normalised'
                        P_dataset.attrs['description'] = 'Parameter covariance'

                        res_dataset = k_group.create_dataset('res', data=res)
                        res_dataset.attrs['unit'] = 'volt'
                        res_dataset.attrs['description'] = 'Radiation field residual'

                        E_theta_dataset = k_group.create_dataset('E_theta', data=E_est_theta)
                        E_theta_dataset.attrs['unit'] = 'volt'
                        E_theta_dataset.attrs['description'] = 'Estimated radiation field polar component'

                        E_phi_dataset = k_group.create_dataset('E_phi', data=E_est_phi)
                        E_phi_dataset.attrs['unit'] = 'volt'
                        E_phi_dataset.attrs['description'] = 'Estimated radiation field azimuthal component'

                        line = f'{tkp:5d}\t{_k:2d}\t{residual_energy:.3e}\t{norm_trace_cov:.3e}'

                        print(line)

        # ========== ========== ========== ========== ========== ==========




        return

        q, P, res, max_rel_error, \
            theta, phi, E_theta, E_phi, \
            E_est_theta, E_est_phi = \
            rfd.estimate_mode_coeff_unbiased_and_regularised(10, 100)

        print(f'mean residual: {res.mean()}')
        print(max_rel_error, P.diagonal().max())

        rfd.plot_mode_power_distribution(q)