import binascii
import glob, os


def check_start_end(hex_line):

    first_chars = hex_line[0:8]
    last_chars1 = hex_line[-8::]
    cal = hex_line[-2::]
    calm = hex_line[-4:-2]
    cala = hex_line[-6:-4]
    calg = hex_line[-8:-6]
    if (first_chars == '0d0a6400') and (last_chars1.isdigit()):
        if int(calm) + int(cala) + int(calg) < 10:
            return True
        else:
            return False
    else:
        return False


def HexView(filein, fileout):
    with open(filein,
              'rb') as in_file , open(fileout, "wb") \
               as outF:

        while True:
            hexdata = in_file.read(100).hex() # I like to read 100 bytes in then new line it.

            if check_start_end(hexdata) is True:
                outF.write(binascii.unhexlify(hexdata))
            else:
                h = in_file.read(1).hex()

                while h != '0d':
                    h = in_file.read(1).hex()
                    # print(h)
                    temp = []
                    if h == '0d':
                        temp = h
                        hd = in_file.read(99).hex()
                        temp += hd

                        if check_start_end(temp) is True:
                            hexdata = temp
                            # print(hexdata)
                            outF.write(binascii.unhexlify(hexdata))
                            break
                        else:
                            h = in_file.read(1)
                    if len(h) == 0:  # breaks loop once no more binary data is read
                        break
            if len(hexdata) == 0:  # breaks loop once no more binary data is read
                break

        # content = in_file.read()
        # print(hexdata)


def run():
    os.chdir('/Users/georgecowie/Documents/Master/Masteroppgave/Master_thesis/MatLab/B13/')
    for file in glob.glob("*.txt"):
        outfile = 'r_' +file
        HexView(file,outfile)


if __name__ == '__main__':
    run()
