import pyvisa
import time
import csv
import os


def run_seed_power_sweep():
    rm = pyvisa.ResourceManager()

    # --- 1. FOLDER CONFIGURATION ---
    new_folder = "Power_Sweep_Test_01"
    os.makedirs(new_folder, exist_ok=True)

    santec_ip = "192.168.1.100"
    osa_ip = "192.168.55.10"

    # --- 2. SEED LASER PARAMETERS ---
    fixed_seed_wl = 1545.0  # Fix the Seed wavelength

    start_power = 0.0  # Initial power in dBm
    end_power = 20.0  # Final power in dBm
    step_size = 1.0  # Power step size in dBm

    # --- 3. OSA PARAMETERS ---
    center_wl = 1549.646
    osa_span = 90.0
    osa_res = 0.1

    try:
        print("Starting sequence for the Power Sweep...")

        # --- SANTEC SEED LASER CONFIG ---
        print("Connecting to Santec Seed Laser...")
        seed = rm.open_resource(f"TCPIP0::{santec_ip}::5000::SOCKET")
        seed.write_termination = "\r"
        seed.read_termination = "\r"
        seed.timeout = 5000

        seed.write("*CLS")
        seed.write(":POW:STAT 1")
        time.sleep(0.5)

        # Configure the fixed wavelength before entering the loop
        print(f"Setting Seed wavelength to {fixed_seed_wl} nm...")
        seed.write(f":WAV {fixed_seed_wl:.3f}")
        seed.query("*OPC?")

        # --- YOKOGAWA OSA CONFIG ---
        print("Connecting to Yokogawa OSA...")
        osa = rm.open_resource(f"TCPIP0::{osa_ip}::10001::SOCKET")
        osa.write_termination = "\r\n"
        osa.read_termination = "\n"
        osa.timeout = 60000

        osa.query('open "anonymous"')
        osa.query(" ")

        print(f"Configuring OSA: Center {center_wl} nm, Span {osa_span} nm...")
        osa.write(f":SENS:WAV:CENT {center_wl}nm")
        osa.write(f":SENS:WAV:SPAN {osa_span}nm")
        osa.write(f":SENS:BWID:RES {osa_res}nm")

        print(f"\nStarting FWM Power Sweep. Data will be saved in: ./{new_folder}/")
        current_power = start_power

        while round(current_power, 2) <= round(end_power, 2):
            pow_str = f"{current_power:.2f}"
            print(f"\n---STEP: Seed Power at {pow_str} dBm ---")

            # 1. Modify Seed power
            seed.write(f":POW {pow_str}")
            seed.query("*OPC?")

            # Extra time to stabilize the laser cavity and internal mechanical attenuator
            time.sleep(1.5)

            # 2. Force Trace A to accept new data
            osa.write(":TRAC:ATTR:TRA WRIT")

            # 3. Trigger the OSA
            print("Acquiring OSA trace...")
            osa.write(":INIT:SMOD SING")
            osa.write(":INIT")

            # 4. THE GOLDEN DELAY
            time.sleep(0.5)

            # 5. Wait for sweep completion
            osa.query("*OPC?")

            # 6. Download data
            print("Downloading data arrays...")
            raw_x = osa.query(":TRAC:DATA:X? TRA")
            raw_y = osa.query(":TRAC:DATA:Y? TRA")

            wls = [float(val) for val in raw_x.split(",")]
            pws = [float(val) for val in raw_y.split(",")]

            # 7. Save to CSV (Modified to reflect power in filename)
            filename = f"fwm_power_sweep_{pow_str}dBm.csv"
            filepath = os.path.join(new_folder, filename)

            with open(filepath, mode="w", newline="") as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(["Wavelength [m]", "Power [dBm]"])
                for x, y in zip(wls, pws):
                    writer.writerow([x, y])

            print(f"Saved {len(wls)} points to: {filepath}")

            # Increment power for the next step
            current_power += step_size

        print("\nPower Sweep completed successfully.")

    except pyvisa.errors.VisaIOError as e:
        print(f"\nWARNING: Network or communication error: {e}")
    except Exception as e:
        print(f"\nWARNING: Unexpected system error: {e}")

    finally:
        print("\nClosing connections...")
        if "seed" in locals():
            seed.close()
        if "osa" in locals():
            try:
                osa.write("close")
                time.sleep(0.2)
            except:
                pass
            osa.close()


if __name__ == "__main__":
    run_seed_power_sweep()
