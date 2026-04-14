import pyvisa
import time
import csv
import os

def run_attenuator_sweep():
    rm = pyvisa.ResourceManager()

    # --- 1. CONFIGURACIÓ DE CARPETES ---
    new_folder = "Attenuator_Sweep_Test_03"
    os.makedirs(new_folder, exist_ok=True)

    santec_ip = "192.168.1.100"
    osa_ip = "192.168.55.10"

    # --- 2. PARÀMETRES DEL SEED LASER (ATENUADOR) ---
    fixed_seed_wl = 1545.0  # Fixem la longitud d'ona del Seed
    
    # Valors de l'atenuador en dB. (Nota: Més atenuació = Menys potència final)
    # Ex: Si start_att = 3.0 (aprox 10.4 dBm) i end_att = 13.0 (aprox 0.4 dBm)
    start_att = 1.0       
    end_att = 15.0        
    step_size = 1.0       

    # --- 3. PARÀMETRES DE L'OSA ---
    center_wl = 1549.646
    osa_span = 27.0
    osa_res = 0.1

    try:
        print("Iniciant seqüència per a l'Attenuator Sweep...")

        # --- SANTEC SEED LASER CONFIG ---
        print("Connectant al Santec Seed Laser...")
        seed = rm.open_resource(f"TCPIP0::{santec_ip}::5000::SOCKET")
        seed.write_termination = "\r"
        seed.read_termination = "\r"
        seed.timeout = 5000

        seed.write("*CLS")
        
        # Desactivem l'Auto Power Control (ALC) per assegurar que fa cas a l'atenuador
        try:
            seed.write(":POW:ALC 0") 
        except:
            pass # Si el model no té aquesta comanda, l'ignorem

        seed.write(":POW:STAT 1")
        time.sleep(0.5)

        # Configurem la longitud d'ona fixa abans d'entrar al bucle
        print(f"Fixant longitud d'ona del Seed a {fixed_seed_wl} nm...")
        seed.write(f":WAV {fixed_seed_wl:.3f}")
        seed.query("*OPC?")

        # --- YOKOGAWA OSA CONFIG ---
        print("Connectant al Yokogawa OSA...")
        osa = rm.open_resource(f"TCPIP0::{osa_ip}::10001::SOCKET")
        osa.write_termination = "\r\n"
        osa.read_termination = "\n"
        osa.timeout = 60000

        osa.query('open "anonymous"')
        osa.query(" ")

        print(f"Configurant OSA: Center {center_wl} nm, Span {osa_span} nm...")
        osa.write(f":SENS:WAV:CENT {center_wl}nm")
        osa.write(f":SENS:WAV:SPAN {osa_span}nm")
        osa.write(f":SENS:BWID:RES {osa_res}nm")

        print(f"\nIniciant FWM Sweep. Les dades es desaran a: ./{new_folder}/")
        current_att = start_att

        while round(current_att, 2) <= round(end_att, 2):
            att_str = f"{current_att:.2f}"
            print(f"\n---STEP: Atenuador a {att_str} dB ---")

            # 1. Modificar el nivell d'atenuació
            seed.write(f":POW:ATT {att_str}")
            seed.query("*OPC?")

            # Temps extra perquè el motor mecànic de l'atenuador acabi de moure's
            time.sleep(1.5)

            # 2. Forçar la Traça A a acceptar noves dades
            osa.write(":TRAC:ATTR:TRA WRIT")

            # 3. Disparar l'OSA
            print("Adquirint traça de l'OSA...")
            osa.write(":INIT:SMOD SING")
            osa.write(":INIT")

            # 4. EL RETARD D'OR (THE GOLDEN DELAY)
            time.sleep(0.5)

            # 5. Esperar que acabi l'escombrada
            osa.query("*OPC?")

            # 6. Descarregar dades
            print("Descarregant matrius de dades...")
            raw_x = osa.query(":TRAC:DATA:X? TRA")
            raw_y = osa.query(":TRAC:DATA:Y? TRA")

            wls = [float(val) for val in raw_x.split(",")]
            pws = [float(val) for val in raw_y.split(",")]

            # 7. Desar en CSV (El nom reflecteix l'atenuació, no la potència final)
            filename = f"fwm_att_sweep_{att_str}dB.csv"
            filepath = os.path.join(new_folder, filename)

            with open(filepath, mode="w", newline="") as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow(["Wavelength [m]", "Power [dBm]"])
                for x, y in zip(wls, pws):
                    writer.writerow([x, y])

            print(f"S'han desat {len(wls)} punts a: {filepath}")

            # Incrementar l'atenuació pel següent pas (això baixarà la potència de sortida)
            current_att += step_size

        print("\nAttenuator Sweep completat amb èxit.")

    except pyvisa.errors.VisaIOError as e:
        print(f"\nALERTA: Error de xarxa o comunicació: {e}")
    except Exception as e:
        print(f"\nALERTA: Error inesperat del sistema: {e}")

    finally:
        print("\nTancant connexions...")
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
    run_attenuator_sweep()
