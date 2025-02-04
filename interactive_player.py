import Interactive_Module as IM

k = IM.query_num("\nEnter the number of colours: \n>", "Enter a single number")
n = IM.query_num("\nEnter the length of the sequence: \n>", "Enter a single number")

hidden = IM.generate(k,n)

print("\n\nLet's start the game!")
tries = 1
while tries<=10:
    round = IM.play(k,n,hidden)
    if round[0] == len(hidden):
        print(f"\nNice! You have guessed the correct sequence in {tries} tries!")
        break
    print(f"You guessed {round[0]} positions correctly \namong the missed ones, {round[1]} of the colours you gave are on a wrong place \nTake another try :)")
    
    tries += 1
    


