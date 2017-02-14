from nose.tools import assert_equal, assert_true

from flash.flmake.utils import hash_list_to_dict, hash_dict_to_str

def test_hash_list_to_dict1():
    obs = hash_list_to_dict(range(4))
    exp = {0: {1: {2: {3: {}}}}}
    assert_equal(obs, exp)


def test_hash_list_to_dict2():
    obs = hash_list_to_dict(range(4))
    obs = hash_list_to_dict(range(0, 8, 2), obs)
    exp = {0: {1: {2: {3: {}}}, 2: {4: {6: {}}}}}
    assert_equal(obs, exp)


def test_hash_dict_to_str1():
    hd = {'0': {'1': {'2': {'3': {}}}}}
    obs = hash_dict_to_str(hd)
    exp = ("+-0\n"
           "  +-1\n"
           "    +-2\n"
           "      +-3\n"
           )
    assert_equal(obs, exp)


def test_hash_dict_to_str2():
    hd = {'0': {'1': {'2': {'3': {}}}, '2': {'4': {'6': {}}}}}
    obs = hash_dict_to_str(hd)
    exp = ("+-0\n"
           "  +-1\n"
           "  | +-2\n"
           "  |   +-3\n"
           "  +-2\n"
           "    +-4\n"
           "      +-6\n"
           )
    assert_equal(obs, exp)


def test_hash_dict_to_str3():
    hd = {'0': {'1': {'2': {'3': {}}}}}
    labels = {'0': 'zero', '2': 'two'}
    obs = hash_dict_to_str(hd, labels)
    exp = ("+-zero (0)\n"
           "  +-1\n"
           "    +-two (2)\n"
           "      +-3\n"
           )
    assert_equal(obs, exp)


