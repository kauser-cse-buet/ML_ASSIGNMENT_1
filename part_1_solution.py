import random

# get number of one present in the binary number v.
def one(v):
    return bin(v).count("1")

# get_random_binary_string_n_bit returns a random binary value of n bits.
def get_random_binary_string_n_bit(n):
    return random.randint(0, 2 ** n - 1)

# toggle_bit_at_position returns a binary number same as bin_number but toggled the bit at position n. Such as: bin_number = 0011,
# toggle_bit_at_position(0011, 0) will return 0010.
def toggle_bit_at_position(bin_number, n):
    musk = (1 << n)
    return bin_number ^ musk

# provided function for the assignment 1, part 1.
def f(v):
    return abs(12 * one(v) - 160)

# implemented algorithm of iterated hill climber.
def iterated_hillclimber(maximum, bit_length_v):
    for t in range(maximum):
        local = False
        v_c = get_random_binary_string_n_bit(bit_length_v)

        while True:
            v_n = 0

            for n in range(bit_length_v):
                v = toggle_bit_at_position(bin_number=v_c, n=n)

                if n == 0:
                    v_n = v
                    # for n = 0, as there are no v  in past so, v_n will be set as v.
                else:
                    if f(v) > f(v_n):
                        # we have to select the v for which f is maximum. that's why we set v_n = v when we found one.
                        v_n = v

            if f(v_c) < f(v_n):
                # if the selected toggled v_n , the maximum value of f for it is greater than that of v_c, that means we get new maximum. so, we update our current maximum.
                v_c = v_n
            else:
                # If no new maximum is found, we have to reset. That's why we set local to True
                local = True

            if local == True:
                print f(v_c),',',
                break

if __name__ == "__main__":
    iterated_hillclimber(100, 40)