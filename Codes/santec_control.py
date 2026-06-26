import pyvisa
import time
import sys

def control_santec_attenuation():
    rm = pyvisa.ResourceManager()
    santec_ip = "192.168.1.100"

    try:
        print("WASSUP GANG")
        print(f"Booting remote attenuation controller for Santec at {santec_ip}...")
        
        # --- SANTEC CONNECTION ---
        seed = rm.open_resource(f"TCPIP0::{santec_ip}::5000::SOCKET")
        seed.write_termination = "\r"  
        seed.read_termination = "\r"   
        seed.timeout = 5000

        # Initialize and ensure laser is ON
        seed.write("*CLS")  
        seed.write(":POW:STAT 1")  
        time.sleep(0.5)
        
        print("\n" + "="*45)
        print("   SANTEC LASER REMOTE CONTROLLER (ATTENUATOR)")
        print("="*45)
        print("Type a number to set the ATTENUATION in dB.")
        print("  (0 dB = maximum power out, higher dB = less power)")
        print("Type 'q' or 'exit' to quit.")
        print("="*45)

        # --- INTERACTIVE CONTROL LOOP ---
        while True:
            user_input = input("\nEnter desired attenuation (dB) [or 'q' to quit]: ").strip().lower()
            
            if user_input in ['q', 'quit', 'exit']:
                print("Ending remote session...")
                break
                
            try:
                # Convert string to float to ensure it's a valid number
                att_val = float(user_input)
                
                # Sanity check (attenuation usually can't be negative)
                if att_val < 0:
                    print("⚠️ Attenuation cannot be negative. Setting to 0 dB.")
                    att_val = 0.0
                
                print(f"Setting attenuation to {att_val:.2f} dB...")
                
                # Command to change the attenuator value
                # Note: If your specific Santec model uses a different SCPI command for the attenuator, 
                # change ":POW:ATT" to that command (e.g., just ":ATT")
                seed.write(f":POW:ATT {att_val:.2f}")
                
                # Ask the laser if it finished executing the command
                seed.query("*OPC?")
                print(f"✅ Attenuation successfully updated to {att_val:.2f} dB!")
                
            except ValueError:
                print("❌ Invalid input! Please enter a numeric value (e.g., 5.0 or 12.5).")

    except pyvisa.errors.VisaIOError as e:
        print(f"\nALERTA MEDUSIL: Network/Communication error: {e}")
    except Exception as e:
        print(f"\nALERTA MEDUSIL: Unexpected system error: {e}")

    finally:
        print("\nClosing connection...")
        if "seed" in locals():
            # Optional: Turn off the laser when exiting
            # seed.write(":POW:STAT 0") 
            seed.close()
            print("Connection closed. See you in the streets.")

if __name__ == "__main__":
    control_santec_attenuation()
