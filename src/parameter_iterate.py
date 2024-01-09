"""
iterate through points in [-1,1]^d
with each coordinate
of the form n/d^k organized by smaller
values of k first (n is not a multiple of d)
when d>1, treat each point by it's maximum k among the coordinates
"""

from typing import Callable, Optional, Tuple, TypeVar, Iterable
# pylint:disable=import-error
from adic_rational import AdicRational
from iter_utils import make_n_param_version

A = TypeVar("A")

class OneParameterIterator:
    """
    one parameter iterator
    """
    def __init__(self, denominator_localization : int,
                 hard_denominator_power_cut : Optional[int] = None):
        assert denominator_localization > 1
        self.denominator_localization = denominator_localization
        self.current_denominator = 1
        self.current_denominator_power = 0
        self.hard_denominator_power_cut = hard_denominator_power_cut
    def __iter__(self):
        """
        itself
        """
        return self
    def __next__(self) -> Tuple[Iterable[int],Tuple[int,int]]:
        """
        all those with the current denominator
        then multiply the denominator by denominator_localization
        to go to the next value of k
        """
        if self.hard_denominator_power_cut is not None and \
            self.current_denominator_power>self.hard_denominator_power_cut:
            raise StopIteration
        if self.current_denominator == 1:
            self.current_denominator_power += 1
            self.current_denominator *= self.denominator_localization
            numerators : Iterable[int] = [-1,0,1]
            denominator_power = 0
        else:
            numerators = (i for i in range(-self.current_denominator,self.current_denominator)\
                          if i % self.denominator_localization != 0)
            denominator_power = self.current_denominator_power
            self.current_denominator_power += 1
            self.current_denominator *= self.denominator_localization
        return (numerators,(self.denominator_localization,denominator_power))

class MultiParameterIterator:
    """
    [-1,1]^d cube parameter iterator
    """
    def __init__(self, denominator_localization : int, num_dimensions : int,
                 expected_per_denominator : Optional[Callable[[int],Optional[int]]] = None,
                 hard_denominator_power_cut : Optional[int] = None):
        assert denominator_localization > 1
        assert num_dimensions >= 1
        self._num_dimensions = num_dimensions
        self._one_param = OneParameterIterator(denominator_localization, hard_denominator_power_cut)
        self._underlying = \
            make_n_param_version(
                map(
                    lambda fracs : \
                        (AdicRational(num,fracs[1][0],fracs[1][1]) for num in fracs[0]),
                    self._one_param
                ),
                self._num_dimensions
            )
        if expected_per_denominator is not None:
            self.__decimate_me(expected_per_denominator)

    def __decimate_me(self,expected_per_denominator : Callable[[int],Optional[int]]):
        def decimator(tup : Tuple[AdicRational,...]) -> bool:
            cur_denom_power = max(map(lambda r: r.denominator_power,tup))
            how_many_should_be = expected_per_denominator(cur_denom_power)
            if how_many_should_be is None:
                return True
            #pylint:disable=unused-variable
            cur_denom_count_leq,cur_denom_count_equal = \
                AdicRational(1,tup[0].denominator_base,cur_denom_power).count_with_denom()
            raise NotImplementedError
            #pylint:disable=unreachable
            how_many_present = 0
            if how_many_should_be >= how_many_present:
                return True
        self._underlying = filter(decimator,self._underlying)

    def __iter__(self):
        """
        itself
        """
        return self

    def __next__(self) -> Tuple[AdicRational,...]:
        """
        the next point as a tuple of num_dimensions floats
        """
        return self._underlying.__next__()

if __name__ == "__main__":
    p = OneParameterIterator(2)
    for cur in p:
        # pylint:disable = invalid-name
        to_print = ""
        for num in cur[0]:
            to_print += f"{num}/{cur[1][0]}^{cur[1][1]},"
        print(to_print)
        if p.current_denominator>16:
            break
    p = OneParameterIterator(4,0)
    np = make_n_param_version(
        map(lambda fracs : (num/(fracs[1][0]**fracs[1][1]) for num in fracs[0]),p),3)
    # pylint:disable = invalid-name
    index = 0
    for cur in np:
        print(cur)
        index += 1
        if index>40:
            break
    else:
        print("Finished the iterator before needing to break manually")
    print("Try again")
    np = MultiParameterIterator(2,2,lambda _: 9)
    # pylint:disable = invalid-name
    index = 0
    for cur in np:
        print(cur)
        index += 1
        if index>80:
            break
