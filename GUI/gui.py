import tkinter as tk
from tkinter import messagebox, scrolledtext
import subprocess
import threading

def run_script():
    # Check if the parameters.dat file exists and execute the shell script
    try:
        with open("parameters.dat", "r") as f:
            params = f.read()
            
        display_terminal_output("Running the script...\n")

        # Run the shell script in a separate thread to keep the UI responsive
        thread = threading.Thread(target=execute_script)
        thread.start()
        
    except FileNotFoundError:
        messagebox.showerror("Error", "parameters.dat file not found!")
    except Exception as e:
        messagebox.showerror("Error", str(e))
        
def execute_script():
    # Function to execute the shell script and capture its output
    process = subprocess.Popen(
        ["bash", "your_script.sh"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    for line in process.stdout:
        display_terminal_output(line)

    for line in process.stderr:
        display_terminal_output(line)

    process.stdout.close()
    process.stderr.close()
    process.wait()

    display_terminal_output("Script finished.\n")

def display_terminal_output(message):
    terminal_output.insert(tk.END, message)
    terminal_output.see(tk.END)  # Scroll to the end of the text box

def show_format():
    format_text = (
        "The 'parameters.dat' file should have the following format:\n\n"
        "xlen = <float>\n"
        "ylen = <float>\n"
        "zlen = <float>\n"
        "rho = <float>\n"
        "uinf = <float>\n"
        "re = <float>\n"
        "ret = <float>\n"
        "tStart = <float>\n"
        "dt = <float>\n"
        "split_t = <integer>\n"
        "zplus = <list of integers>\n"
        "chplot = <string>\n"
        "in_path = <string>\n"
        "out_path = <string>\n"
        "rans_path = <string>\n\n"
        "Example:\n"
        "xlen = 1.0\n"
        "ylen = 2.0\n"
        "zlen = 3.0\n"
        "rho = 1.225\n"
        "uinf = 0.5\n"
        "re = 10000\n"
        "ret = 200\n"
        "tStart = 0.0\n"
        "dt = 0.01\n"
        "split_t = 10\n"
        "zplus = [10, 20, 30]\n"
        "chplot = 'normal' or 'all'\n"
        "in_path = /path/to/input_directory\n"
        "out_path = /path/to/output_directory\n"
        "rans_path = /path/to/rans_datas_file"
    )
    messagebox.showinfo("Parameter File Format", format_text)
    
def show_parameters():
    # Display the contents of parameters.dat file
    try:
        with open("parameters.dat", "r") as f:
            params = f.read()
        messagebox.showinfo("Current Parameters", params)
    except FileNotFoundError:
        messagebox.showerror("Error", "parameters.dat file not found!")
    except Exception as e:
        messagebox.showerror("Error", str(e))

# Create the main window
root = tk.Tk()
root.title("Run Simulation")

# Styling variables
font = ("Arial", 14)
bg_color = "#e6f2ff"  # Light blue background
button_color = "#4da6ff"  # Blue button
text_color = "#003366"  # Dark blue text

# Set the window background color
root.configure(bg=bg_color)

# Create and place labels and buttons
tk.Label(root, text="Run the cflow analysis with the 'parameters.dat' file", font=font, bg=bg_color, fg=text_color).pack(pady=20)

button_format = tk.Button(root, text="Show File Format", font=font, command=show_format)
button_format.pack(pady=10)

button_show_params = tk.Button(root, text="Show Parameters", font=font, command=show_parameters)
button_show_params.pack(pady=10)

button_run = tk.Button(root, text="Run", font=font, bg=button_color, fg="white", command=run_script)
button_run.pack(pady=10)

# Add a ScrolledText widget to display terminal output
terminal_output = scrolledtext.ScrolledText(root, height=10, font=("Courier", 12), bg="#f2f2f2", fg="black")
terminal_output.pack(pady=20, padx=10, fill=tk.BOTH, expand=True)



# Run the application
root.mainloop()
