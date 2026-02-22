from judge import get_feedback
import random



def generate(k,n):
       
    """
    generates a random sequence for a player to guess

    Args:
        k(int): number of colors considered
        n(int): length of the sequence

    Returns:
        list: sequence of n numbers from 1 to k to be guessed

    """
    hidden = [random.randint(1,k) for x in range(1,n+1)]
    return hidden


def play(k, n, hidden):
    
    """
    Asks a player to enter a sequence, then gives a feedback about correctness

    Args:
        k(int): number of colors considered
        n(int): length of the sequence
        hidden(list): a sequence to be guessed by a player 

    Returns:
        tuple: feedback 

    """
    while True:
        try:
            query = query_seq("\nEnter your guess: \n>", "\n Enter only single numbers separated with single space")
        
        except EOFError:
            print(f"\nIt's a pity you gave up. The hidden sequence was {hidden}")
            exit()
        if correct_format(k, hidden, query):
            feedback = get_feedback(query, hidden)
            return feedback
        else:
            print(f"\nYou need to enter a sequence of {n} numbers from 1 to {k} separated with space")
            return play(k,n,hidden)
        

    


def query_seq(inp, exc):
    
    """
    Asks a player to input a sequence of numbers

    Args:
        inp(str): Message for a player before he enters the input
        exc(str): Message for a player after he enters a wrong input

    Returns:
        list: sequence entered by a player

    Raises: 
        ValueError: {exc}

    """
    while True:
        try:
            answer = input(inp)
            answer = list(map(int, answer.split()))
            return answer
        except ValueError:
            print(exc)


def query_num(inp, exc):
    
    """
    Asks a player to input a single number

    Args:
        inp(str): Message for a player before he enters the input
        exc(str): Message for a player after he enters a wrong input

    Returns:
        list: number entered by a player

    Raises: 
        ValueError: {exc}

    """
    while True:
        try:
            answer = int(input(inp))
            return answer
        except ValueError:
            print(exc)
    


def correct_format(k, hidden, query):
    """
    Checks if a player's input is matching requirements

    Args:
        k(int): Number of colors considered
        hidden(list): A sequence to be guessed by a player
        query(list): The player's query

    Returns:
        bool: True if query format is correct

    """
    correct = [marble for marble in query if marble in range(1,k+1)]
    if len(correct) == len(hidden): 
        return True
    


