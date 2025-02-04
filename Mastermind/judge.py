def get_feedback(guess: list[int], code: list[int]) -> tuple[int,int]:
       
    """
    Calculates how many positions in sequence are guessed correctly and how many colors (values) are correct but on a wrong place.

    Args:
        guess (list): sequence to be compared with code
        code (list): sequence that guess is being compared to

    Returns:
        tuple: (number of positions correctly guessed, number of colours in wrongplace)

    """
    
    correct_place = sum(g == c for g, c in zip(guess, code))
    wrong_place = sum(min(guess.count(c), code.count(c)) for c in set(guess)) - correct_place #ile razy cyferka z query miałaby możliwość być na tym samym miejscu co cyferka w code minus to kiedy rzeczywiście tam jest
    return (correct_place, wrong_place)


    
