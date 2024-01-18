"""
rational of the form a/b^d
"""

from typing import Tuple

class AdicRational:
    """
    unreduced rational number
    where denominator is power of given base
    """
    def __init__(self,numerator : int,denominator_base : int,denominator_power : int):
        self.numerator = numerator
        self.denominator_base = denominator_base
        self.denominator_power = denominator_power

    def __repr__(self):
        if self.denominator_power == 1:
            return f"{self.numerator}/{self.denominator_base}"
        if self.denominator_power == 0:
            return f"{self.numerator}"
        return f"{self.numerator}/{self.denominator_base}^{self.denominator_power}"

    def count_with_denom(self) -> Tuple[int,int]:
        """
        how many have this same denominator power
        how many have denominator powers less than or equal to this one's
        """
        if self.denominator_power == 0:
            return 3,3
        def count_leq(denom_power) -> int:
            top_power = self.denominator_base**denom_power
            return top_power*2 + 1
        def count_equal(denom_power) -> int:
            return self.denominator_base**(denom_power-1)*(self.denominator_base-1)*2
        my_count_leq = count_leq(self.denominator_power)
        my_count_equal = count_equal(self.denominator_power)
        return my_count_leq, my_count_equal

    def to_float(self) -> float:
        """
        convert to float
        """
        return self.numerator/(self.denominator_base**self.denominator_power)
