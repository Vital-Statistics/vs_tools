#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 16:03:23 2024

@author: rudy
"""

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

