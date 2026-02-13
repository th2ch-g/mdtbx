import argparse
import math
from ..logger import generate_logger

LOGGER = generate_logger(__name__)

# Constants
A0 = -59.2194
A1 = 0.07594
B0 = -22.8396
B1 = 0.01347
D0 = 1.1677
D1 = 0.002976
KB = 0.008314  # kJ/mol/K

MAXITER = 100

def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "gen_temperatures",
        help="Generate temperatures for REMD simulations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--pdes", type=float, required=True, help="Exchange probability")
    parser.add_argument("--tlow", type=float, required=True, help="Lower temperature limit")
    parser.add_argument("--thigh", type=float, required=True, help="Upper temperature limit")
    parser.add_argument("--nw", type=int, required=True, help="Number of water molecules")
    parser.add_argument("--np", type=int, required=True, help="Number of protein atoms")
    parser.add_argument("--tol", type=float, default=1e-4, help="Tolerance")

    parser.add_argument(
        "--pc",
        type=int,
        default=0,
        choices=[0, 1, 2],
        help="Constraints in the protein: 0=Fully Flexible, 1=Bonds to hydrogens only, 2=All bonds",
    )

    parser.add_argument(
        "--wc",
        type=int,
        default=0,
        choices=[0, 2, 3],
        help="Constraints in water: 0=Fully Flexible, 2=Flexible angle, 3=Rigid",
    )

    parser.add_argument(
        "--hff",
        type=int,
        default=0,
        choices=[0, 1],
        help="Hydrogens in protein: 0=All H, 1=Polar H",
    )

    parser.add_argument(
        "--vs",
        type=int,
        default=0,
        choices=[0, 1],
        help="Virtual sites in protein: 0=None, 1=Virtual Hydrogen",
    )

    parser.add_argument(
        "--alg",
        type=int,
        default=0,
        choices=[0, 1],
        help="Simulation type: 0=NPT, 1=NVT (Note: only NPT is fully supported)",
    )

def calc_mu(nw, np_val, temp, fener):
    return ((A0 + A1 * temp) * nw + (B0 + B1 * temp) * np_val - temp * fener)

def myeval(m12, s12, cc, u):
    argument = -cc * u - (u - m12) ** 2 / (2 * s12 * s12)
    return math.exp(argument)

def myintegral(m12, s12, cc):
    int_val = 0.0
    umax = m12 + 5 * s12
    du = umax / 100.0

    u = 0.0
    while u < umax:
        di = myeval(m12, s12, cc, u + du / 2.0)
        int_val += di
        u += du

    return du * int_val / (s12 * math.sqrt(2 * math.pi))

def run(args):
    # Input validation
    if not (0 <= args.pdes <= 1):
        raise ValueError("Exchange probability Pdes must be between 0 and 1")
    if args.thigh <= args.tlow:
        raise ValueError("Upper temperature limit must be greater than lower limit")
    if args.tlow <= 0 or args.thigh <= 0:
        raise ValueError("Temperatures must be > 0")
    if args.np == 0:
        raise ValueError("Number of protein atoms cannot be zero")
    if args.alg != 0:
        # LOGGER.warning("Only NPT (Alg=0) is fully supported as per original implementation.")
        raise ValueError("Can not do constant volume yet!")

    # Calculate derived variables
    npp = 0
    nprot = 0
    nh = 0
    vc = 0
    nc = 0

    if args.hff == 0:
        nh = round(args.np * 0.5134)
        if args.vs == 1:
            vc = round(1.91 * nh)
        nprot = args.np
    else:
        npp = round(args.np / 0.65957)
        nh = round(args.np * 0.22)
        if args.vs == 1:
            vc = round(args.np + 1.91 * nh)
        nprot = npp

    if args.pc == 1:
        nc = nh
    elif args.pc == 2:
        nc = args.np

    ndf = (9 - args.wc) * args.nw + 3 * args.np - nc - vc
    flex_ener = 0.5 * KB * (nc + vc + args.wc * args.nw)

    print("Summary of input and derived variables:")
    print(f"Pdes: {args.pdes}")
    print(f"Temperature range: {args.tlow} - {args.thigh}")
    print(f"Number of water molecules: {args.nw}")
    print(f"Number of protein atoms: {args.np}")
    if npp > 0:
        print(f"Including all H: ~ {npp}")
    print(f"Number of hydrogens in protein: ~ {nh}")
    print(f"Number of constraints: ~ {nc}")
    print(f"Number of vsites: ~ {vc}")
    print(f"Number of degrees of freedom: ~ {ndf}")
    print(f"Energy loss due to constraints: {flex_ener:.2f} (kJ/mol K)")
    print("-" * 40)

    # Main loop
    t_list = [args.tlow]
    p_list = []
    mu_list = []
    sigma_list = []
    mm_list = []
    ss_list = []

    # Initial values for T1 (first temperature)
    mu_list.append(calc_mu(args.nw, nprot, args.tlow, flex_ener))
    sigma_list.append(math.sqrt(ndf) * (D0 + D1 * args.tlow))

    while t_list[-1] < args.thigh:
        t1 = t_list[-1]
        t2 = t1 + 1.0
        if t2 >= args.thigh:
            t2 = args.thigh

        low = t1
        high = args.thigh

        iter_count = 0
        piter = 0.0
        forward = True

        while (abs(args.pdes - piter) > args.tol) and (iter_count < MAXITER):
            iter_count += 1
            mu12 = (t2 - t1) * ((A1 * args.nw) + (B1 * nprot) - flex_ener)

            cc = (1.0 / KB) * ((1.0 / t1) - (1.0 / t2))

            var = ndf * (D1 * D1 * (t1 * t1 + t2 * t2) + 2 * D1 * D0 * (t1 + t2) + 2 * D0 * D0)
            sig12 = math.sqrt(var)

            if sig12 == 0:
                raise ValueError("Sigma = 0")

            # I1
            erfarg1 = mu12 / (sig12 * math.sqrt(2))
            i1 = 0.5 * math.erfc(erfarg1)

            # I2
            i2 = myintegral(mu12, sig12, cc)
            piter = i1 + i2

            if piter > args.pdes:
                if forward:
                    t2 = t2 + 1.0
                else:
                    low = t2
                    t2 = low + ((high - low) / 2.0)

                if t2 >= args.thigh:
                    t2 = args.thigh
            elif piter < args.pdes:
                if forward:
                    forward = False
                    low = t2 - 1.0
                high = t2
                t2 = low + ((high - low) / 2.0)

        # Store results for this interval
        p_list.append(piter)
        mm_list.append(mu12)
        ss_list.append(sig12)

        t_list.append(t2)

        # Values for new temperature
        mu_list.append(calc_mu(args.nw, nprot, t2, flex_ener))
        sigma_list.append(math.sqrt(ndf) * (D0 + D1 * t2))

    # Output table
    print("\nTemperatures and Energies")
    print(f"{'Index':<5} {'Temp (K)':<10} {'mu':<10} {'sigma':<10} {'mu12':<10} {'sigma12':<10} {'P12':<10}")

    for k in range(len(t_list)):
        idx = k + 1
        temp = t_list[k]
        mu = mu_list[k]
        sigma = sigma_list[k]

        if k == 0:
             print(f"{idx:<5} {temp:<10.2f} {mu:<10.0f} {sigma:<10.2f}")
        else:
             # Previous interval values
             prev_mu12 = mm_list[k-1]
             prev_sig12 = ss_list[k-1]
             prev_p = p_list[k-1]
             print(f"{idx:<5} {temp:<10.2f} {mu:<10.0f} {sigma:<10.2f} {prev_mu12:<10.1f} {prev_sig12:<10.2f} {prev_p:<10.4f}")

    print("\nTemperature list for scripting:")
    print(", ".join([f"{t:.2f}" for t in t_list]))
