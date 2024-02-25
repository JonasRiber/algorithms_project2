def choose_gapcost():
    while True:
        lin_or_aff = input("Linear/affine gapcost: ").lower()
        if lin_or_aff == "linear":
            user_input = input("Gap cost: ")
            if user_input.isdigit():
                gap_cost = 5 * int(user_input)
                break
            print("Invalid input. Please try again.")
        elif lin_or_aff == "affine":
            user_input = input("Gap cost: ")
            if user_input.isdigit():
                gap_cost = 5 + 5 * int(user_input)
                break
            print("Invalid input. Please try again.")
        print("Invalid input. Please try again.")