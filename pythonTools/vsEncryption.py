"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""
# Auto-generated from vsEncryption/ package modules.

# ---- source: vsEncryption/encrypt_file.py ----
#!/usr/bin/env python3


# Function to encrypt data
def encrypt_file(file_path, key=None):
    import os
    import base64
    from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
    from cryptography.hazmat.primitives import padding
    from cryptography.hazmat.backends import default_backend
    
    if key is None:
        key=base64.b64decode(os.environ['ENCRYPTION_KEY'])
    with open(file_path, 'rb') as file:
        file_data = file.read()
    iv = os.urandom(16)
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend())
    encryptor = cipher.encryptor()
    padder = padding.PKCS7(algorithms.AES.block_size).padder()
    padded_data = padder.update(file_data) + padder.finalize()
    encrypted_data = encryptor.update(padded_data) + encryptor.finalize()
    with open(file_path + '.enc', 'wb') as encrypted_file:
        encrypted_file.write(iv + encrypted_data)


# ---- source: vsEncryption/loadEncrypted.py ----
#!/usr/bin/env python3


import pandas as pd


# Function to decrypt data in memory
def loadEncrypted(encrypted_file_path, key=None):
    from io import BytesIO
    import os
    import base64
    from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
    from cryptography.hazmat.primitives import padding
    from cryptography.hazmat.backends import default_backend
    
    if key is None:
        key=base64.b64decode(os.environ['ENCRYPTION_KEY'])
    with open(encrypted_file_path, 'rb') as encrypted_file:
        iv = encrypted_file.read(16)
        encrypted_data = encrypted_file.read()
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend())
    decryptor = cipher.decryptor()
    padded_data = decryptor.update(encrypted_data) + decryptor.finalize()
    unpadder = padding.PKCS7(algorithms.AES.block_size).unpadder()
    file_data = unpadder.update(padded_data) + unpadder.finalize()
    dDat=BytesIO(file_data)
    dDat.seek(0)
    res=pd.read_parquet(dDat)
    return res
