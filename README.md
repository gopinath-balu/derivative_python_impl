# derivative_python_impl
Python implementation for the derivatives


# TVC(T)

```
import math

def TVC(T, M, M_star, A, h, H, P, theta, alpha, Y, SP, IE, IP):
    """
    Calculate Total Variable Cost (TVC) based on the provided mathematical formula.

    Parameters:
    T (float): Current value of T
    M (float): Parameter M
    M_star (float): Parameter M*
    A (float): Parameter A
    h (float): Parameter h
    H (float): Parameter H
    P (float): Parameter P
    theta (float): Parameter theta
    alpha (float): Parameter alpha
    Y (float): Parameter Y
    SP (float): Parameter SP
    IE (float): Parameter IE
    IP (float): Parameter IP

    Returns:
    float: The computed TVC value.
    """
    # Precompute recurring terms to improve efficiency
    alpha_theta = alpha * theta
    alpha_plus_theta = alpha + theta
    theta_alpha_plus_theta = theta * alpha_plus_theta

    # Calculate terms
    term1 = A / (T + h * H)
    term2 = (2 * P / theta) - (2 * P * T / theta)
    term3 = (P * T / theta) - (1 / (alpha * theta))
    term4 = - (T / theta) + (1 / (alpha * alpha_plus_theta))
    term5 = (2 * T / alpha_plus_theta) - (1 / theta_alpha_plus_theta)
    term6 = - (T / alpha_plus_theta) - (alpha * T / theta_alpha_plus_theta)
    term7 = - (alpha * T / alpha_plus_theta) + (alpha / alpha_plus_theta)
    term8 = (T**2 + (2 * h * T) / theta_alpha_plus_theta)
    term9 = (2 * h * alpha / theta_alpha_plus_theta)
    term10 = (2 * h / alpha_plus_theta) - (h * T / (T * alpha_plus_theta))
    term11 = -(h * alpha * T / alpha_plus_theta) - (h * T * theta / alpha_plus_theta)
    term12 = (2 * P / theta) + (2 * P * T / theta)
    term13 = (Y * P / (T * theta)) + Y - (Y * P / (T * theta))
    term14 = - Y * P + (Y * H / T)
    term15 = - (SP * IE * M / (2 * T)) - (SP * IE * M)
    term16 = Y * IP * (P / theta)
    term17 = - (Y / T) * IP * (theta - alpha_theta * M)
    term18 = (2 * alpha_theta * T) - (alpha_theta * M * T) + (alpha_theta * T / theta_alpha_plus_theta)
    term19 = (2 * P / theta) - (2 * P * M / theta)
    term20 = (1 + alpha * M / (alpha * alpha_plus_theta)) - (P * T / theta)

    # Sum all terms
    TVC_result = (term1 + term2 + term3 + term4 + term5 + term6 + term7 +
                  term8 + term9 + term10 + term11 + term12 + term13 +
                  term14 + term15 + term16 + term17 + term18 + term19 +
                  term20)
    
    return TVC_result

# Example usage with mock data
if __name__ == "__main__":
    T = 10
    M = 5
    M_star = 20
    A = 2
    h = 1
    H = 1
    P = 3
    theta = 2
    alpha = 0.5
    Y = 4
    SP = 2
    IE = 1.5
    IP = 0.8

    result = TVC(T, M, M_star, A, h, H, P, theta, alpha, Y, SP, IE, IP)
    print(f"The computed TVC value is: {result:.2f}")

```

The code has been refined for better efficiency and readability. Key improvements include:

1.  **Precomputing Reusable Terms**: Recurring calculations like - \(\alpha \cdot \theta\)
- \(\theta \cdot (\alpha + \theta)\) are computed once and reused.
2.  **Organized Comments**: Added detailed comments for parameters, calculations, and the function's purpose.
3.  **Efficiency**: Reduced redundant operations.

# First order derivative TVC(T)

```
def dTVC_dT(T, M, A, h, H, P, theta, alpha, Y, SP, IE, IP, e):
    """
    Calculate the first derivative of Total Variable Cost (TVC) with respect to T.

    Parameters:
    T (float): Current value of T
    M (float): Parameter M
    A (float): Parameter A
    h (float): Parameter h
    H (float): Parameter H
    P (float): Parameter P
    theta (float): Parameter theta
    alpha (float): Parameter alpha
    Y (float): Parameter Y
    SP (float): Parameter SP
    IE (float): Parameter IE
    IP (float): Parameter IP
    e (float): Parameter e

    Returns:
    float: The computed value of the first derivative of TVC.
    """
    # Precompute recurring terms for efficiency
    alpha_theta = alpha * theta
    alpha_plus_theta = alpha + theta
    theta_alpha_plus_theta = theta * alpha_plus_theta

    # Calculate terms
    term1 = -A / T
    term2 = (P / theta) - (1 / theta) - (1 / alpha_plus_theta)
    term3 = - (alpha / theta_alpha_plus_theta)
    term4 = - (2 * alpha * T / alpha_plus_theta)
    term5 = - (Y * H / T) + (SP * IE * M / (2 * T))
    term6 = (Y * IP / (T * alpha * alpha_plus_theta)) - (Y * IP * alpha / alpha_plus_theta)
    term7 = (P / theta) - (T / (T * h * H))
    term8 = - (h / alpha_plus_theta) + (alpha * T / alpha_plus_theta)
    term9 = (h * H / T) - (h / (T * alpha_plus_theta))
    term10 = 2 / alpha_plus_theta - (2 * P / theta)
    term11 = -(h * alpha / alpha_plus_theta) - (h * theta / alpha_plus_theta)
    term12 = (alpha * T / alpha_plus_theta)

    # Nested term calculation
    nested_term = (theta * e + P * H * theta + alpha / alpha_plus_theta - e * theta * P + 
                   alpha / theta_alpha_plus_theta + alpha * T / alpha_plus_theta - H * theta + 
                   H * theta * alpha * T - P * T * alpha / (theta * (H - e)) -
                   2 * theta * (H - e) * alpha / alpha_plus_theta + alpha / alpha_plus_theta)

    # Combine all terms
    derivative = (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + 
                  term9 + term10 + term11 + term12) * nested_term

    return derivative

# Example usage with mock data
if __name__ == "__main__":
    T = 10
    M = 5
    A = 2
    h = 1
    H = 1
    P = 3
    theta = 2
    alpha = 0.5
    Y = 4
    SP = 2
    IE = 1.5
    IP = 0.8
    e = 0.1

    result = dTVC_dT(T, M, A, h, H, P, theta, alpha, Y, SP, IE, IP, e)
    print(f"The computed first derivative of TVC is: {result:.2f}")

```

The formula for the first derivative of TVC(T) has been converted into Python code. Key features include:

1.  **Precomputed Terms**: Recurring terms like - \(\alpha \cdot \theta\)
- \(\alpha + \theta\) are calculated once for efficiency.
2.  **Nested Term Handling**: The nested terms are encapsulated into a separate calculation for clarity.
3.  **Parameterization**: The function accepts all parameters as inputs, making it flexible for various scenarios.


# Third order derivative TVC(T)

```
def d2TVC_dT2(T, M, A, h, H, P, theta, alpha, Y, SP, IE, IP, e):
    """
    Calculate the second derivative of Total Variable Cost (TVC) with respect to T.

    Parameters:
    T (float): Current value of T
    M (float): Parameter M
    A (float): Parameter A
    h (float): Parameter h
    H (float): Parameter H
    P (float): Parameter P
    theta (float): Parameter theta
    alpha (float): Parameter alpha
    Y (float): Parameter Y
    SP (float): Parameter SP
    IE (float): Parameter IE
    IP (float): Parameter IP
    e (float): Parameter e

    Returns:
    float: The computed value of the second derivative of TVC.
    """
    # Precompute recurring terms for efficiency
    alpha_plus_theta = alpha + theta
    theta_alpha_plus_theta = theta * alpha_plus_theta
    alpha_theta = alpha * theta

    # Calculate individual terms
    term1 = A / T
    term2 = (2 * h / theta_alpha_plus_theta) + Y * H
    term3 = -(Y * IP / (alpha * alpha_plus_theta)) + (Y * IP * M / alpha_plus_theta)
    term4 = (alpha / alpha_plus_theta) - (h / alpha_plus_theta)

    nested_term1 = (theta * e + P * H * theta + alpha / alpha_plus_theta - e * theta * P +
                    alpha / theta_alpha_plus_theta + alpha * T / alpha_plus_theta - H * theta +
                    H * theta * alpha * T - P * T * alpha / (theta * (H - e)) -
                    2 * theta * (H - e) * alpha / alpha_plus_theta + alpha / alpha_plus_theta)

    term5 = nested_term1

    term6 = (2 * T / (T * h * H)) - (h / alpha_plus_theta)
    term7 = -(2 * alpha / alpha_plus_theta) - (SP * IE * M / (2 * T))
    term8 = (h * H / T) - (h / (T * alpha_plus_theta)) + (2 / alpha_plus_theta) - (2 * P / theta)
    term9 = -(h * alpha / alpha_plus_theta) - (h * theta / alpha_plus_theta)
    term10 = (alpha * T / alpha_plus_theta)

    # Nested term calculation
    nested_term2 = ((theta * (H - e) - 2 * theta * (H - e) * alpha / alpha_plus_theta +
                     2 * theta * e - alpha * e * theta - theta * alpha * e * P + alpha * e /
                     alpha_plus_theta - theta * alpha / alpha_plus_theta * e * (1 + alpha * T) +
                     H * theta * alpha * e - H * theta * alpha * e * (1 + alpha * T) +
                     P * alpha * theta * e * (1 + alpha * T)))

    term11 = nested_term2

    # Combine all terms
    second_derivative = (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 +
                         term9 + term10 + term11)

    return second_derivative

# Example usage with mock data
if __name__ == "__main__":
    T = 10
    M = 5
    A = 2
    h = 1
    H = 1
    P = 3
    theta = 2
    alpha = 0.5
    Y = 4
    SP = 2
    IE = 1.5
    IP = 0.8
    e = 0.1

    result = d2TVC_dT2(T, M, A, h, H, P, theta, alpha, Y, SP, IE, IP, e)
    print(f"The computed second derivative of TVC is: {result:.2f}")

```

The second derivative of TVC(T) has been implemented in Python. Key features include:

1.  **Efficient Precomputations**: Common terms are computed once and reused to save computational resources.
2.  **Clarity in Nested Terms**: Separate calculations for nested terms improve readability and reduce errors.
3.  **Generalization**: All parameters are treated as inputs for flexibility.
