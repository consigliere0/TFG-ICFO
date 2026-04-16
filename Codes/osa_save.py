import pyvisa
import csv
import os
import sys

def acquire_single_trace(pump_power_dbm):
    rm = pyvisa.ResourceManager()
    osa_ip = "192.168.55.10"
    
    # Canviem el nom de la carpeta perquè no es barregin els experiments
    folder = "Pump_Sweep_02"
    os.makedirs(folder, exist_ok=True)

    try:
        print(f"Connectant a l'OSA per descarregar la dada a {pump_power_dbm} dBm...")
        osa = rm.open_resource(f"TCPIP0::{osa_ip}::10001::SOCKET")
        osa.write_termination = "\r\n"
        osa.read_termination = "\n"
        osa.timeout = 10000

        osa.query('open "anonymous"')
        osa.query(" ")

        # Descarregar directament el que hi hagi a la pantalla
        print("Descarregant dades de la traça A...")
        raw_x = osa.query(":TRAC:DATA:X? TRA")
        raw_y = osa.query(":TRAC:DATA:Y? TRA")

        wls = [float(val) for val in raw_x.split(",")]
        pws = [float(val) for val in raw_y.split(",")]

        # Guardar amb el nom de la potència del pump
        filename = f"pump_power_sweep_{pump_power_dbm}dBm.csv"
        filepath = os.path.join(folder, filename)

        with open(filepath, mode="w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Wavelength [m]", "Power [dBm]"])
            for x, y in zip(wls, pws):
                writer.writerow([x, y])
                
        print(f"[OK] Dades guardades a {filepath}")

    except Exception as e:
        print(f"Error: {e}")
    finally:
        if 'osa' in locals():
            osa.close()

if __name__ == "__main__":
    # Comprovem que se li passa un argument per terminal
    if len(sys.argv) < 2:
        print("Ús: python guarda_osa_power.py <Potencia_Pump_dBm>")
        print("Exemple: python guarda_osa_power.py 16.0")
    else:
        power = sys.argv[1]
        acquire_single_trace(power)
