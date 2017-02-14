"""By importing this module, common derivied fields are registered with yt."""

from yt.data_objects.field_info_container import add_field

def _nion(field, data):
    return data['dens'] * data['sumy'] * 6.022E23

add_field ('nion', function=_nion, take_log=True)


def _abar(field, data):
    return 1.0 / data['sumy']

add_field ('abar', function=_abar, take_log=False)


def _velo(field, data):
    return (data['velx']**2 + data['vely']**2 + data['velz']**2)**0.5

add_field ('velo', function=_velo, take_log=True)

add_field ('velr_rz', function=lambda field, data: data['velx'], take_log=True)
add_field ('velz_rz', function=lambda field, data: data['vely'], take_log=True)
add_field ('velt_rz', function=lambda field, data: data['velz'], take_log=True)
