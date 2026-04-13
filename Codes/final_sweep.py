import pyvisa
import time
import csv
import os

def run_fwm_sweep():
    rm = pyvisa.ResourceManager()

    new_folder = "FWM_sweep_test_01"
    os.makedirs(new_folder, exist_ok=True)  

    santec_ip = "192.168.1.100"
    osa_ip = "192.168.55.10"

    # Seed parameters
    start_wl = 1540.0
    end_wl = 1570.0
    step_size = 1.0

    # OSA parameters
    center_wl = 1549.646
    osa_span = 90.0
    osa_res = 0.1

    try:
        print("WASSUP GANG")
        print("Time to conquer these streets. Booting sequence...")
        
        # --- SANTEC SEED LASER CONFIG ---
        print("Connecting to Santec Seed Laser...")
        seed = rm.open_resource(f"TCPIP0::{santec_ip}::5000::SOCKET")
        seed.write_termination = "\r"  
        seed.read_termination = "\r"   
        seed.timeout = 5000

        seed.write("*CLS")  
        seed.write(":POW:STAT 1")  
        time.sleep(0.5)
        seed.write(":POW 10.00")  
        seed.query("*OPC?")  

        # --- YOKOGAWA OSA CONFIG ---
        print("Connecting to Yokogawa OSA...")
        osa = rm.open_resource(f"TCPIP0::{osa_ip}::10001::SOCKET")
        osa.write_termination = "\r\n"  
        osa.read_termination = "\n"
        
        # Ensure the PyVISA timeout is long enough for a full sweep (60 seconds)
        osa.timeout = 60000 

        osa.query('open "anonymous"')  
        osa.query(" ")

        print(f"Configuring OSA: center {center_wl} nm, Span {osa_span} nm...")
        osa.write(f":SENS:WAV:CENT {center_wl}nm")
        osa.write(f":SENS:WAV:SPAN {osa_span}nm")
        osa.write(f":SENS:BWID:RES {osa_res}nm")

        print(f"\nStarting FWM sweep. Data will be saved in: ./{new_folder}/")
        current_wl = start_wl

        while round(current_wl, 3) <= round(end_wl, 3):
            wl_str = f"{current_wl:.3f}"
            print(f"\n---STEP: Seed at {wl_str} nm ---")

            # 1. Move seed
            seed.write(f":WAV {wl_str}")
            seed.query("*OPC?")
            time.sleep(1.5) 

            # 2. Force Trace A to accept new data
            osa.write(":TRAC:ATTR:TRA WRIT")

            # 3. Trigger OSA
            print("Acquiring OSA trace...")
            osa.write(":INIT:SMOD SING")
            osa.write(":INIT")

            # 4. THE GOLDEN DELAY
            # Give the OSA's CPU time to physically engage the sweep motor
            time.sleep(0.5)

            # 5. Wait for sweep completion
            # Now that the motor is running, *OPC? will hold Python here until it finishes
            osa.query("*OPC?")

            # 6. Download data
            print("Downloading array data...")
            raw_x = osa.query(":TRAC:DATA:X? TRA")
            raw_y = osa.query(":TRAC:DATA:Y? TRA")

            wls = [float(val) for val in raw_x.split(",")]
            pws = [float(val) for val in raw_y.split(",")]

            # 7. Save to CSV
            filename = f"fwm_sweep_{wl_str}.csv"
            filepath = os.path.join(new_folder, filename)

            with open(filepath, mode="w", newline="") as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(["Wavelength [m]", "Power [dBm]"])
                for x, y in zip(wls, pws):
                    writer.writerow([x, y])
                    
            print(f"Saved {len(wls)} points to: {filepath}")
            current_wl += step_size

        print("\nWavelength sweep completed WIIIII.")

    except pyvisa.errors.VisaIOError as e:
        print(f"\nALERTA MEDUSIL: Network/Communication error: {e}")
    except Exception as e:
        print(f"\nALERTA MEDUSIL: Unexpected system error: {e}")

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
    run_fwm_sweep()
